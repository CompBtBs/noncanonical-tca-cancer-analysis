#load the libraries
import pandas as pd
import re
import numpy as np
from progressbar import ProgressBar,Bar,Percentage
from scanpy import AnnData
from cobra.flux_analysis.variability import find_essential_reactions,find_essential_genes

import os
from lxml import etree as ET
from colour import Color
from IPython.display import Image
import seaborn as sns
from PIL import Image as PilImage
from  matplotlib.colors import to_hex
from pdf2image import convert_from_path
from PIL import  ImageFont, ImageDraw 
"""
Class to compute the RAS values

"""

class RAS_computation:

    def __init__(self,adata,model):
                                                       
        self._logic_operators = ['and', 'or', '(', ')']
        self.val_nan = np.nan

        # Build the dictionary for the GPRs
        df_reactions = pd.DataFrame(index=[reaction.id for reaction in model.reactions])
        gene_rules=[reaction.gene_reaction_rule for reaction in model.reactions]   
        
        gene_rules=[el.replace("OR","or").replace("AND","and").replace("(","( ").replace(")"," )") for el in gene_rules]        
        df_reactions['rule'] = gene_rules
        df_reactions = df_reactions.reset_index()
        df_reactions = df_reactions.groupby('rule').agg(lambda x: sorted(list(x)))
        
        self.dict_rule_reactions = df_reactions.to_dict()['index']

        # build useful structures for RAS computation
        self.model = model
        self.count_adata = adata.copy()
        self.genes = self.count_adata.var.index.intersection([gene.id for gene in model.genes])
        
        #check if there is one gene at least 
        if len(self.genes)==0:
            print("ERROR: no gene of the count matrix is in the metabolic model!")
            print(" are you sure that the gene annotation is the same for the model and the count matrix?")
            return -1
        
        self.cell_ids = list(self.count_adata.obs.index.values)
        self.count_df_filtered = self.count_adata.to_df().T.loc[self.genes]
 
    def compute(self,
                or_expression=np.nansum,    # type of operation to do in case of an or expression (max, sum, mean)
                and_expression=np.nanmin,   # type of operation to do in case of an and expression(min, sum)
                drop_na_rows=True,          # if True remove the nan rows of the ras  matrix
                drop_duplicates=False,      # if true, remove duplicates rows
                regexp=re.compile(r"\([a-zA-Z0-9-.:\s]+\)"),  # regular expression inside a parenthesis
                print_progressbar=True,     # if True, print the progress bar
                add_count_metadata=True,    # if True add metadata of cells in the ras adata
                add_met_metadata=True,      # if True add metadata from the metabolic model (gpr and compartments of reactions)
                add_essential_reactions=False,
                add_essential_genes=False
                ):

        self.or_function = or_expression
        self.and_function = and_expression
        
        ras_df = pd.DataFrame(index=range(len(self.dict_rule_reactions)), columns=self.cell_ids)
        ras_df[:][:] = self.val_nan
        
        if print_progressbar:
            pbar = ProgressBar(widgets=[Percentage(), Bar()], maxval=len(self.dict_rule_reactions)).start()
            i = 0
        
        # for loop on reactions
        ind = 0       
        for rule, reaction_ids in self.dict_rule_reactions.items():
            if len(rule) != 0:
                # there is one gene at least in the formula
                rule_split = rule.split()
                rule_split_elements = list(filter(lambda x: x not in self._logic_operators, rule_split))  # remove of all logical operators
                rule_split_elements = list(np.unique(rule_split_elements))                                # genes in formula
                
                # which genes are in the count matrix?                
                genes_in_count_matrix = list(set([el for el in rule_split_elements if el in self.genes]))
                genes_notin_count_matrix = list(set([el for el in rule_split_elements if el not in self.genes]))


                if len(genes_in_count_matrix) > 0: #there is at least one gene in the count matrix
                     if len(rule_split) == 1:
                         #one gene --> one reaction
                         ras_df.iloc[ind] = self.count_df_filtered.loc[genes_in_count_matrix]
                     else:                        
                        # more genes in the formula
                        lista = re.findall(regexp, rule)
                        if len(lista) == 0:
                             #or/and sequence
                             matrix = self.count_df_filtered.loc[genes_in_count_matrix].values
                             if len(genes_notin_count_matrix) > 0:
                                matrix = np.vstack([matrix, [self.val_nan for el in self.cell_ids]])

                             if 'or' in rule_split: 
                                ras_df.iloc[ind] = self.or_function(matrix, axis=0)
                             else:
                                ras_df.iloc[ind] = self.and_function(matrix, axis=0)
                        else:
                            # ho almeno una tonda
                            data = self.count_df_filtered.loc[genes_in_count_matrix]  # dataframe of genes in the GPRs
                            genes = data.index
                            j = 0
                             
                            for cellid in self.cell_ids:    #for loop on the cells
                                lista_cell = lista.copy()
                                rule_cell = rule
                                 
                                while len(lista_cell) > 0:
                                    #
                                    for el in lista_cell:
                                        #print(el[1:-1])
                                        value = self._evaluate_expression(el[1:-1].split(), data[cellid], genes)
                                        rule_cell = rule_cell.replace(el, str(value))   
                                    lista_cell = re.findall(regexp, rule_cell)      
         
                                ras_df.iloc[ind, j] = self._evaluate_expression(rule_cell.split(), data[cellid], genes)
                                j=j+1
      
            ind = ind+1
            #update percentage
            if print_progressbar:
                pbar.update(i+1)
                i = i+1
        
        if print_progressbar:
            pbar.finish()
        
        ras_df=ras_df.astype("float")    
        ras_df['REACTIONS'] = [reaction_ids for rule,reaction_ids in self.dict_rule_reactions.items()]
        
        reactions_common = pd.DataFrame()
        reactions_common["REACTIONS"] = ras_df['REACTIONS']
        reactions_common["proof2"] = ras_df['REACTIONS']
        reactions_common = reactions_common.explode('REACTIONS')
        reactions_common = reactions_common.set_index("REACTIONS")

        ras_df = ras_df.explode("REACTIONS")
        ras_df = ras_df.set_index("REACTIONS")

        if drop_na_rows:
            ras_df = ras_df.dropna(how="all")
            
        if drop_duplicates:
            ras_df = ras_df.drop_duplicates()
        
        #create AnnData structure for RAS
        ras_adata = AnnData(ras_df.T)

        #add metadata
        if add_count_metadata:
            ras_adata.var["common_gprs"] = reactions_common.loc[ras_df.index]
            ras_adata.var["common_gprs"] = ras_adata.var["common_gprs"].apply(lambda x: ",".join(x))
            for el in self.count_adata.obs.columns:
                ras_adata.obs["countmatrix_"+el]=self.count_adata.obs[el]

        if add_met_metadata:
            if len(self.model.compartments)>0:
                  ras_adata.var['compartments']=[list(self.model.reactions.get_by_id(reaction).compartments) for reaction in ras_adata.var.index]  
                  ras_adata.var['compartments']=ras_adata.var["compartments"].apply(lambda x: ",".join(x))
            
            ras_adata.var['GPR rule'] = [self.model.reactions.get_by_id(reaction).gene_reaction_rule for reaction in ras_adata.var.index]

        if add_essential_reactions:            
            essential_reactions=find_essential_reactions(self.model)
            essential_reactions=[el.id for el in essential_reactions]            
            ras_adata.var['essential reactions']=["yes" if el in essential_reactions else "no" for el in ras_adata.var.index]
        
        if add_essential_genes:
            essential_genes=find_essential_genes(self.model)
            essential_genes=[el.id for el in essential_genes]
            ras_adata.var['essential genes']=[" ".join([gene for gene in genes.split()  if gene in essential_genes]) for genes in ras_adata.var["GPR rule"]]
        
        return ras_adata


    def _check_number(self,value):
      try:
        float(value)
        return True
      except ValueError:
        return False

    def _evaluate_expression(self, rule_split, values_cellid, genes):
        
        #ci sono per forza solo or
        rule_split2 = list(filter(lambda x: x != "or" and x!="and", rule_split))   

        values = list()
        i=0
        for el in rule_split2:
             if self._check_number(el):
                 values.append(float(el))
             elif el in genes:
                 values.append(values_cellid[el])
             else:
                 values.append(self.val_nan)
                 i=i+1
                 
        if i==len(rule_split2):
            return self.val_nan
        if "or" in rule_split:
            #or sequence
            return self.or_function(values)
        else:
            #and sequence
            return self.and_function(values)



class map_network():

    def __init__(self):
        return None
    
    def fix_style(self, col, width):
        tmp = []
        tmp.append('stroke:' + col)
        tmp.append('stroke-width:' + str(float(width)))
        tmp.append('stroke-opacity:1')
        tmp.append('fill-opacity:1')
        return ';'.join(tmp)
        
    def colorMapRAS(self,
                     fileMap,
                     fileMapColor,
                     df,                            #dataset of reaction to color the map
                     colormap=["blue", "red"],        #color scale              #
                     nosignificant_color="grey",    #color of reaction whose differences results not significant
                     ref_width_value=10,            #dimension of reaction rows
                     nosignificant_width=5,         #dimension of the row of reaction whose differences results not significant
                     min_val=0,                     #minimum width value. The minimum width of the row is min_val*ref_width_value
                     max_val=10,                    #maximum width value. The maximum width of the row is max_val*ref_width_value
                     unconfined=True,               #Set unconfined=True to disable max-width confinement of the image.
                     width_image=1000               #dimension of the image
                     ):

            if not isinstance(df, pd.DataFrame):
                df=pd.DataFrame(df)

            df=df[df[df.columns[0]].isna() == False]
            
            #imposto quelle da colorare
            df_significant = df[np.abs(df[df.columns[0]]) > 0]
            df_significant = df_significant.sort_values(by=df_significant.columns[0])
            
            #imposto quelle non significative
            df_nosignificant = df[np.abs(df[df.columns[0]]) == 0]

            #create a vector
            color_range = [Color(colormap[0]),Color(colormap[1])]

            dict_colors = dict()
            dict_widths = dict()
            
            for el in df_significant.index:
                dict_colors[el] = str(color_range[0]) if df_significant.loc[el,df.columns[0]]<0 else str(color_range[1])
                
                width=np.max([np.abs(df_significant.loc[el,df.columns[0]]),min_val])
                width=np.min([width,max_val])
                
                dict_widths[el] = str(np.round(ref_width_value*width,2))
                    

            root = ET.parse(open(fileMap, 'r'))


            for key in dict_colors.keys():
                for element in root.iter():          
                    if element.get('id') == "R_" + key:
                        stringa = element.get("style")
                        if stringa is not None:
                            # Modifica lo stile esistente
                            stringa = stringa.replace("#000000", dict_colors[key])
                            stringa2 = stringa.split(";")
                            stringa3 = list()
                            for el in stringa2:
                                if "stroke-width:" in el:
                                    stringa3.append("stroke-width:" + str(dict_widths[key]))
                                else:
                                    stringa3.append(el)
                            element.set("style", ";".join(stringa3))
                        else:
                            element.set("style", self.fix_style(dict_colors[key], dict_widths[key]))

            with open(fileMapColor, 'w') as newMapFile:
                string=ET.tostring(root, pretty_print=True)
                string=self._to_utf8(string)
                newMapFile.write(string)

            if width_image:
                image = Image(url=fileMapColor, unconfined =unconfined, width=width_image)
            else:
                image = Image(url=fileMapColor, unconfined =unconfined)

            return image

   
    def colorMapUpDown(self,
                 fileMap,
                 fileMapColor,
                 dfTotal,                       #dataset of reaction to color the map
                 colormap=["blue","red"],
                 absolute_value=True,
                 nosignificant_color="grey",    #color of reaction whose differences results not significant
                 ref_width_value=10,            #dimension of reaction rows
                 nosignificant_width=5,         #dimension of the row of reaction whose differences results not significant
                 max_val=None,                    #maximum width value. The maximum width of the row is max_val*ref_width_value
                 unconfined=True,               #Set unconfined=True to disable max-width confinement of the image.
                 width_image=1000,              #dimension of the image
                 center_range=False,
                 scale_factor=2,
                 thres_binary=0,
                 title_text=None,
                 dimension=(6005,15),
                 title_font=None
                 ):

        #dfTotal deve avere score e diagnostica e la voce reverse_test e reverse
        df=pd.DataFrame(dfTotal["score"])
        df=df[df[df.columns[0]].isna() == False]
        
         
        #quali sono significative (con dentro anche quelle di diagnostica)
        df_significant = df[np.abs(df[df.columns[0]]) > 0]
        df_significant = df_significant.sort_values(by=df_significant.columns[0])

        #quali reazioni non sono significative
        df_nosignificant = df[np.abs(df[df.columns[0]]) == 0]

        dict_colors = dict()
        dict_widths = dict()       

        min_value=df_significant.min()
        max_value=df_significant.max() 

        if max_val==None:
            if (np.abs(min_value)>np.abs(max_value)).values:
                max_val=np.abs(min_value.values)
            else:
                max_val=np.abs(max_value.values)
        
        df_significant_colors=df_significant.copy()
        df_significant_width=df_significant.copy()
        df_significant_width=df_significant_width.abs()

        if center_range:
            df_significant_colors=np.round((df_significant_colors-max_value)/max_val*100)
            df_significant_width=np.round(df_significant_colors/max_val,4)
        else:
            df_significant_colors=np.round((df_significant_colors)/max_value*100)
            df_significant_width=np.round(df_significant_width/max_val*scale_factor,4)

        #setto il colore e lo spessore delle significative
        for el in df_significant_colors.index:
            #valore del colore
            
            #print(index)
            if dfTotal.loc[el,"reverse_test"]==False:
                #entrambe positive
                if df_significant_colors.loc[el,df.columns[0]]<thres_binary:
                    dict_colors[el]=str(Color(colormap[0])) 
                else:
                    dict_colors[el]=str(Color(colormap[1]))
                    
            elif dfTotal.loc[el,"reverse_test"]==True:
                #entrambe negative
                if df_significant_colors.loc[el,df.columns[0]]>thres_binary:
                    dict_colors[el]=str(Color(colormap[0])) 
                else:
                    dict_colors[el]=str(Color(colormap[1]))                    #entrambe negative
                
            elif dfTotal.loc[el,"reverse_test"]=="+to-":
                dict_colors[el]="#00FFFF"  
                dict_widths[el]=str(nosignificant_width)
                
            elif dfTotal.loc[el,"reverse_test"]=="-to+":
                dict_colors[el]="green" #d+
                dict_widths[el]=str(nosignificant_width)
            #valore larghezza
            width=df_significant_width.loc[el,df.columns[0]]
            dict_widths[el] = str(np.round(ref_width_value*width,4))
                
        #setto il colore e lo spessore delle non significative
        for el in df_nosignificant.index:
            dict_colors[el]=nosignificant_color
            dict_widths[el]=str(nosignificant_width)
            
        root = ET.parse(open(fileMap, 'r'))

        #qui si colora
        for key in dict_colors.keys():
            for element in root.iter():
                if element.get('id') == "R_" + key:
                    stringa = element.get("style")
                    if stringa is not None:
                        # Modifica lo stile esistente
                        stringa = stringa.replace("#000000", dict_colors[key])
                        stringa2 = stringa.split(";")
                        stringa3 = list()
                        for el in stringa2:
                            if "stroke-width:" in el:
                                stringa3.append("stroke-width:" + str(dict_widths[key]))
                            else:
                                stringa3.append(el)
                        element.set("style", ";".join(stringa3))
                    else:
                        element.set("style", self.fix_style(dict_colors[key], dict_widths[key]))


                if  element.get('id')=="F_"+key and dfTotal.loc[key,"reverse_test"]==False:
                    stringa = element.get("style")
                    #caso standard positivo e positivo F
                    if stringa is not None:
                        stringa3=self.replace_string(element,dict_colors[key],dict_widths[key])
                        element.set("style",";".join(stringa3))
                    else:
                        element.set("style", self.fix_style(dict_colors[key], dict_widths[key]))

                if  element.get('id')=="B_"+key and dfTotal.loc[key,"reverse_test"]==False:
                    stringa = element.get("style")
                    if stringa is not None: 
                        #caso standard positivo e positivo B
                        stringa3=self.delete_string(element,dict_colors[key],dict_widths[key])
                        element.set("style",";".join(stringa3))
                    else:
                        element.set("style", self.fix_style("#00000", "0"))

                if  element.get('id')=="B_"+key and dfTotal.loc[key,"reverse_test"]==True:
                    stringa = element.get("style") 
                    if stringa is not None:
                        #ribalto la freccia
                        stringa3=self.replace_string(element,dict_colors[key],dict_widths[key])
                        element.set("style",";".join(stringa3))
                    else:
                        element.set("style", self.fix_style(dict_colors[key], dict_widths[key]))

                if  element.get('id')=="F_"+key and dfTotal.loc[key,"reverse_test"]==True: 
                    stringa = element.get("style")
                    if stringa is not None:
                        #caso standard positivo e positivo B
                        stringa3=self.delete_string(element,dict_colors[key],dict_widths[key])
                        element.set("style",";".join(stringa3))
                    else:
                        element.set("style", self.fix_style("#00000", "0"))
                if  element.get('id')=="F_"+key and (dfTotal.loc[key,"reverse_test"]=="+to-" or dfTotal.loc[key,"reverse_test"]=="-to+"): 
                    stringa = element.get("style")
                    if stringa is not None:
                        #ribalto la freccia
                        stringa3=self.replace_string(element,dict_colors[key],dict_widths[key])
                        element.set("style",";".join(stringa3))
                    else:
                        element.set("style", self.fix_style(dict_colors[key], dict_widths[key]))

                      
        with open(fileMapColor, 'w') as newMapFile:
            string=ET.tostring(root, pretty_print=True)
            string=self._to_utf8(string)
            newMapFile.write(string)

        if width_image:
            image = Image(url=fileMapColor, unconfined =unconfined, width=width_image)
        else:
            image = Image(url=fileMapColor, unconfined =unconfined)
        

        
        return image,dict_colors,dict_widths

    def colorMapCentroid(self,
                     fileMap,
                     fileMapColor,
                     df,                            #dataset of reaction to color the map
                     colormap="magma",        #color scale #["blue","red"]
                     reverse_colormap=True,
                     nosignificant_color="grey",    #color of reaction whose differences results not significant
                     ref_width_value=10,            #dimension of reaction rows
                     nosignificant_width=2,         #dimension of the row of reaction whose differences results not significant
                     max_val=None,                    #maximum width value. The maximum width of the row is max_val*ref_width_value
                     unconfined=True,               #Set unconfined=True to disable max-width confinement of the image.
                     width_image=1000,              #dimension of the image
                     minimum_width=1,
                     max_value=None,
                     min_value=None,
                     apply_log=True,
                     scale_factor=2,
                     thres_binary=0,
                     range_color=5,
                     title_text="",
                     dimension=(6005,15),
                     title_font=None
                 ):

        
        #df è il dataframe originale con segno
        #df_significant è in modulo
        
        if not isinstance(df, pd.DataFrame):
            df=pd.DataFrame(df)
    
        #quali sono significative
        df_significant = df.sort_values(by=df.columns[0])#ordino dal più piccolo al più grande
        df_significant = df_significant[np.abs(df_significant[df_significant.columns[0]]) != 0]
        
        eps=0
        #quali reazioni non sono significative (ovvero 0)
        df_nosignificant = df[np.abs(df[df.columns[0]])<= eps]

        #apply log scale
        if apply_log:
            if max_value:
                max_value=np.log1p(eps+np.abs(max_value))
            if min_value:
                min_value=np.log1p(eps+np.abs(min_value))
            df_significant[df_significant.columns[0]]=df_significant[df_significant.columns[0]].apply(lambda x: np.log1p(eps+np.abs(x)))

        else:
            df_significant[df_significant.columns[0]]=df_significant[df_significant.columns[0]].apply(lambda x: eps+np.abs(x))         
        
        
        df_significant=df_significant[df_significant[df_significant.columns[0]].isna() == False]

        dict_colors = dict()
        dict_widths = dict()       
        
        if type(colormap)==str:
            try:
                palette = sns.color_palette(colormap, range_color)
            except:
                palette = sns.mpl_palette(colormap, range_color)
            colors=[to_hex(el) for el in palette]
        else:
            colors=[el for el in Color(colormap[0]).range_to(colormap[1],range_color)]             

        if reverse_colormap:
            colors=sorted(colors, reverse=True)        

        if max_value is None:
            max_value=df_significant.abs().max()

        if min_value is None:
            min_value=df_significant.abs().min()
           
        
        range_val=max_value-min_value
        
        df_significant_colors=df_significant.copy()
        df_significant_width=df_significant.copy().abs()
        
        df_significant_colors=np.round((df_significant_colors-min_value)/range_val*(range_color-1))
        df_significant_width=np.round(df_significant_width/max_value*scale_factor,4)
        
        #setto il colore e lo spesso delle significative
        for el in df_significant_colors.index:
            
            index=int(df_significant_colors.loc[el,df.columns[0]])
            dict_colors[el]=str(colors[index])
            #valore larghezza
            width=df_significant_width.loc[el,df.columns[0]]
            dict_widths[el] = str(np.round(ref_width_value*width,4)+minimum_width)
                
        #setto il colore e lo spessore delle non significative
        for el in df_nosignificant.index:
            dict_colors[el]=nosignificant_color
            dict_widths[el]=str(nosignificant_width)
        

        f=open(fileMap, "r")
        contents = f.readlines()
        f.close()
        
        #contents.insert(18, create_rectangle(colors))
        #contents ="".join(contents)
        
        f=open(fileMapColor, "w")
        f.writelines(contents)
        f.close()
        
        
        root= ET.parse(open(fileMapColor, 'r')) 

        for key in dict_colors.keys():
            for element in root.iter():
                #REAZIONE
                if element.get('id') == "R_" + key:
                    stringa = element.get("style")
                    if stringa is not None:
                        # Modifica lo stile esistente
                        stringa = stringa.replace("#000000", dict_colors[key])
                        stringa2 = stringa.split(";")
                        stringa3 = list()
                        for el in stringa2:
                            if "stroke-width:" in el:
                                stringa3.append("stroke-width:" + str(dict_widths[key]))
                            else:
                                stringa3.append(el)
                        element.set("style", ";".join(stringa3))
                    else:
                        element.set("style", self.fix_style(dict_colors[key], dict_widths[key]))
        
                if  element.get('id')=="F_"+key and df.loc[key,df.columns[0]]>0: 
                    stringa = element.get("style")
                    if stringa is not None:
                        # Modifica lo stile esistente
                        stringa = stringa.replace("#000000", dict_colors[key])
                        stringa2 = stringa.split(";")
                        stringa3 = list()
                        for el in stringa2:
                            if "stroke-width:" in el:
                                stringa3.append("stroke-width:" + str(dict_widths[key]))
                            else:
                                stringa3.append(el)
                        element.set("style", ";".join(stringa3))
                    else:
                        element.set("style", self.fix_style(dict_colors[key], dict_widths[key]))

                if  element.get('id')=="B_"+key and df.loc[key,df.columns[0]]>0: 
                    
                    stringa=element.get("style")
                    if stringa is not None:
                        stringa=stringa.replace("stroke-opacity:1","stroke-opacity:0")
                        stringa=stringa.replace("fill-opacity:1","fill-opacity:0")
                        
                        
                        stringa2=stringa.split(";")
                        stringa3=list()
                        for el in stringa2:
                            if "stroke-width:" in el:
                                #ref_width_value=float(el.split(":")[1])
                                
                                stringa3.append("stroke-width:"+str(dict_widths[key]))
                            else:
                                stringa3.append(el)
                            
                        element.set("style",";".join(stringa3))

                    else:
                        # Crea una nuova stringa di stile
                        tmp = self.fix_style(dict_colors[key], dict_widths[key])

                        tmp=tmp.replace("stroke-opacity:1","stroke-opacity:0")
                        tmp=tmp.replace("fill-opacity:1","fill-opacity:0")
                        
                        stringa2=tmp.split(";")
                        stringa3=list()
                        for el in stringa2:
                            if "stroke-width:" in el:
                                #ref_width_value=float(el.split(":")[1])
                                
                                stringa3.append("stroke-width:"+str(dict_widths[key]))
                            else:
                                stringa3.append(el)
                            
                        element.set("style",";".join(stringa3))
                        

                if  element.get('id')=="B_"+key and df.loc[key,df.columns[0]]<0: 

                    stringa = element.get("style")
                    if stringa is not None:
                        # Modifica lo stile esistente
                        stringa = stringa.replace("#000000", dict_colors[key])
                        stringa2 = stringa.split(";")
                        stringa3 = list()
                        for el in stringa2:
                            if "stroke-width:" in el:
                                stringa3.append("stroke-width:" + str(dict_widths[key]))
                            else:
                                stringa3.append(el)
                        element.set("style", ";".join(stringa3))
                    else:
                        element.set("style", self.fix_style(dict_colors[key], dict_widths[key]))


                if  element.get('id')=="F_"+key and df.loc[key,df.columns[0]]<0: 
                    stringa=element.get("style")
                    if stringa is not None:
                        stringa=stringa.replace("stroke-opacity:1","stroke-opacity:0")
                        stringa=stringa.replace("fill-opacity:1","fill-opacity:0")
                        
                        
                        stringa2=stringa.split(";")
                        stringa3=list()
                        for el in stringa2:
                            if "stroke-width:" in el:
                                #ref_width_value=float(el.split(":")[1])
                                
                                stringa3.append("stroke-width:"+str(dict_widths[key]))
                            else:
                                stringa3.append(el)
                            
                        element.set("style",";".join(stringa3))

                    else:
                        # Crea una nuova stringa di stile
                        tmp = self.fix_style(dict_colors[key], dict_widths[key])

                        tmp=tmp.replace("stroke-opacity:1","stroke-opacity:0")
                        tmp=tmp.replace("fill-opacity:1","fill-opacity:0")
                        
                        stringa2=tmp.split(";")
                        stringa3=list()
                        for el in stringa2:
                            if "stroke-width:" in el:
                                #ref_width_value=float(el.split(":")[1])
                                
                                stringa3.append("stroke-width:"+str(dict_widths[key]))
                            else:
                                stringa3.append(el)
                            
                        element.set("style",";".join(stringa3))


                if  ((element.get('id')=="F_"+key) or (element.get('id')=="B_"+key))  and df.loc[key,df.columns[0]]==0: 

                    stringa = element.get("style")
                    if stringa is not None:
                        # Modifica lo stile esistente
                        stringa = stringa.replace("#000000", dict_colors[key])
                        stringa2 = stringa.split(";")
                        stringa3 = list()
                        for el in stringa2:
                            if "stroke-width:" in el:
                                stringa3.append("stroke-width:" + str(dict_widths[key]))
                            else:
                                stringa3.append(el)
                        element.set("style", ";".join(stringa3))
                    else:
                        element.set("style", self.fix_style(dict_colors[key], dict_widths[key]))



        with open(fileMapColor, 'w') as newMapFile:
            string=ET.tostring(root, pretty_print=True)
            string=self._to_utf8(string)
            newMapFile.write(string)
        
        if width_image:
            image = Image(url=fileMapColor, unconfined =unconfined, width=width_image)
        else:
            image = Image(url=fileMapColor, unconfined =unconfined)

        
        return image,dict_colors

    def delete_string(self,element,dict_colors_value,dict_widths_value):
        stringa=element.get("style")

        if stringa is not None:
            stringa=stringa.replace("stroke-opacity:1","stroke-opacity:0")
            stringa=stringa.replace("fill-opacity:1","fill-opacity:0")
            
            stringa2=stringa.split(";")
            stringa3=list()
            for el in stringa2:
                if "stroke-width:" in el:
                    #ref_width_value=float(el.split(":")[1])
                    
                    stringa3.append("stroke-width:"+str(dict_widths_value))
                else:
                    stringa3.append(el)

        else:
            # Crea una nuova stringa di stile
            tmp = self.fix_style(dict_colors_value, dict_widths_value)

            tmp=tmp.replace("stroke-opacity:1","stroke-opacity:0")
            tmp=tmp.replace("fill-opacity:1","fill-opacity:0")
            
            stringa2=tmp.split(";")
            stringa3=list()
            for el in stringa2:
                if "stroke-width:" in el:
                    #ref_width_value=float(el.split(":")[1])
                    
                    stringa3.append("stroke-width:"+str(dict_widths_value))
                else:
                    stringa3.append(el)

        return stringa3
    
                        

    def replace_string(self,element,dict_colors_value,dict_widths_value):
        stringa=element.get("style")
        if stringa is not None:
            stringa=stringa.replace("#000000",dict_colors_value)
            
            stringa2=stringa.split(";")
            stringa3=list()
            for el in stringa2:
                if "stroke-width:" in el:
                    #ref_width_value=float(el.split(":")[1])
                    
                    stringa3.append("stroke-width:"+str(dict_widths_value))
                else:
                    stringa3.append(el)
        else:
            stringa3 = self.fix_style(dict_colors_value, dict_widths_value)
        return stringa3  

    def _to_utf8(self,s):
        return s if isinstance(s, str) else s.decode('utf-8') 
