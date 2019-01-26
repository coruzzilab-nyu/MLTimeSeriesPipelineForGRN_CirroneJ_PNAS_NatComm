import sys

import os

import pandas as pd

import itertools

import matplotlib.pyplot as plt

import numpy as np

from random import shuffle

import networkx as nx

from sklearn.metrics import *

import matplotlib.pyplot as plt


#plt.switch_backend('agg')

class Validation:

    def __init__(self):
        pass
        #Read gene names and description
        gene_desc_filename = "./Genes/Genesnames_descrpt.txt"
        self.gene_desc_df = pd.read_table(gene_desc_filename, header=None)
        self.gene_desc_df[0] = map(lambda x: x.upper(), self.gene_desc_df[0])
        self.gene_desc_df.set_index(0, inplace=True)


    def run(self, params):
        pvalue_thres, list_name, remove_noExpInProtopl, remove_noJITs, \
        cheat_meth, thres_prec, plot_one_random_ordering, all_random_ordering, \
        filename_genes_to_remove_prot, filename_genes_to_remove_noJIT = params 

        if not os.path.exists("./results"):
            os.makedirs("./results")

        #Validate TARGET
        targetspath="./TFTargets/"

        if not os.path.exists(targetspath):
            os.makedirs(targetspath)
            print "No TFs lists in folder: ", targetspath
            print "Please add TFs lists to: ", targetspath
            sys.exit()

            
        tfs_list_tovalidate = []
        num_tfs_to_validate = 0

        for file in os.listdir(targetspath):
            if file.startswith("AT"):
                parts = file.split("_")
                tf_name = parts[0]
                tfs_list_tovalidate.append(tf_name)  
                num_tfs_to_validate = num_tfs_to_validate + 1
                    

                    

        #1658genes 145 TFs
        geneNames_file = "geneNames1658g.txt"
        TFNames_file = "TFNames145.txt"
        matrix_net_dfg_path = "./newDFG/Root1658g145TF_Quant_df2_1e2/"
        wAvg_matrix_name = "wAvgArabidopsis_n1658g145TFs_root.txt"
        tissue = "roots"
        NgenesNtfs = "genes1658_numtfs145"



        matrix_net_dfg = pd.read_table(matrix_net_dfg_path+wAvg_matrix_name, header=None)
        matrix_net_dfg_2 = matrix_net_dfg.loc[:, (matrix_net_dfg != 0).any(axis=0)]

        genenames = pd.read_table(matrix_net_dfg_path+geneNames_file, header=None)
        tfnames = pd.read_table(matrix_net_dfg_path+TFNames_file, header=None)


        #Update list of TFs to validate based on universe of TFs
        tfs_list_tovalidate = list(set(tfs_list_tovalidate).intersection(set(tfnames[0])))
        num_of_Tfs_tovalid = len(tfs_list_tovalidate)

        if num_of_Tfs_tovalid>0:

            valid_set_name = "_RootsOnly_"+str(pvalue_thres)+"Pval_"+str(num_of_Tfs_tovalid)+"TFs_"+list_name

            print "Num of total TFs to validate: ", num_of_Tfs_tovalid

            tfnames_ = list(tfnames.values.tolist())
            tfs_merged_list = list(itertools.chain.from_iterable(tfnames_))
            matrix_net_dfg_2.columns = tfs_merged_list


            genenames_ = list(genenames.values.tolist())
            genenames_merged_list = list(itertools.chain.from_iterable(genenames_))
            matrix_net_dfg_2.index = genenames_merged_list

            if tissue == "roots" and remove_noExpInProtopl:
                num_targets_tmp = matrix_net_dfg_2.shape[0]
                print "Matrix dfg before removing genes: ", matrix_net_dfg_2.shape
                targets_to_remove_file=filename_genes_to_remove_prot
                genes_to_remove = pd.read_table(targets_to_remove_file, header=None)
                genes_to_remove_in_univ = list(set(genenames_merged_list).intersection(set(genes_to_remove[0])))
                matrix_net_dfg_2.drop(genes_to_remove_in_univ, inplace=True)
                print "Matrix Shape after filtering genes not expressed in protoplast: ", matrix_net_dfg_2.shape
                print "Total genes removed: ", str( num_targets_tmp - matrix_net_dfg_2.shape[0])



            if tissue == "roots" and remove_noJITs:
                num_targets_tmp = matrix_net_dfg_2.shape[0]
                print "Matrix dfg before removing genes: ", matrix_net_dfg_2.shape
                targets_to_remove_file=filename_genes_to_remove_noJIT
                genes_to_remove = pd.read_table(targets_to_remove_file, header=None)
                genes_to_remove_in_univ = list(set(matrix_net_dfg_2.index.values).intersection(set(genes_to_remove[0])))
                matrix_net_dfg_2.drop(genes_to_remove_in_univ, inplace=True)
                print "Matrix Shape after filtering genes not expressed in protoplast: ", matrix_net_dfg_2.shape
                print "Total genes removed: ", str( num_targets_tmp - matrix_net_dfg_2.shape[0])


            dfg_edges = matrix_net_dfg_2.transpose().stack().reset_index()
            dfg_edges.columns = ['TF','Target','weight']
            dfg_edges.sort_values(["weight"], ascending=False, inplace=True)  
            dfg_edges.reset_index(inplace=True)
            dfg_edges.drop("index", axis=1, inplace=True)


            dfg_edges.columns = ["source", "dest", "score"]

            dfg_edges.to_csv("./results/FullNetwork_WithSIGN_NewDFGEdges__"+NgenesNtfs+"_datatypeTS_"+tissue+".txt", sep="\t", index=False)


            dfg_edges_abs = dfg_edges.copy()

            dfg_edges_abs['score'] = abs(dfg_edges_abs['score'])

            dfg_edges_abs.sort_values('score', inplace=True, ascending=False)

            dfg_edges_abs.reset_index(inplace=True)
            dfg_edges_abs.drop('index', axis=1, inplace=True)


            #print "####NO SiGN####"
            dfg_edges_abs.to_csv("./results/FullNetwork_NoSIGN_NewDFGEdges__"+NgenesNtfs+"_datatypeTS_"+tissue+".txt", sep="\t", index=False)

            validated_edges_df = self.validateEdges(dfg_edges_abs, tfs_list_tovalidate, pvalue_thres)


            prcurve_namefile = "./results/_PRcurve_"+str(cheat_meth)+"cheatMeth_NewDFGEdges__"+NgenesNtfs+"_datatypeTS_"+tissue+"_GoldStand"+str(num_of_Tfs_tovalid)+str(valid_set_name)+"_nosign"
            network_namefile = "./results/_Network_"+str(cheat_meth)+"cheatMeth_NewDFGEdges__"+NgenesNtfs+"_datatypeTS_"+tissue+"_GoldStand"+str(num_of_Tfs_tovalid)+str(valid_set_name)+"_nosign"


            validated_edges_df_keys, ranking4TFs_filtered0rows0cols, edges_gs, edges_score, precision, recall, average_precision = self.print_PR_curve(tfs_list_tovalidate, dfg_edges_abs, prcurve_namefile, validated_edges_df, cheat_meth, network_namefile, thres_prec, all_random_ordering, plot_one_random_ordering)
            print "AUPR DFG: ", average_precision

        else:
            print "Error No edges in the gold standard"

    def getTargets_fromNetwork_dataframe(self, df, tfname):
        df_tmp = df.set_index("source")
        return df_tmp.loc[tfname]["dest"]

    def interesect_twogenelists(self, genelist1, genelist2):
        return set(genelist1).intersection(set(genelist2))

    def getgene_pvaluefdr(self, df, genename):
        df_tmp = df.set_index("target")
        return str(round(df_tmp.loc[genename]["foldchange"],3))+"/"+str(round(df_tmp.loc[genename]["fdr"],3))


    def getgene_rankingscore(self, df, tfname, target):
        df_tmp = df.set_index(["source", "dest"])
        return round(float(1)*float(df_tmp.loc[tfname,target]["score"]), 1)

    def get_genedescr(self, genename):
        return self.gene_desc_df.loc[genename]



    def validateEdges(self, ranking, tfs_list_tovalidate, pvalue_thres):

        tot_num_edges = ranking.shape[0]
        num_edges_filtered = ranking.shape[0]
        percent_of_edges_filter = 100*(float(num_edges_filtered) / float(tot_num_edges))
        
        # print "Total number of edges before filtering: ", tot_num_edges
        # print "Total number of edges after filtering: ", num_edges_filtered

        # print "Edges' percentage of the total network after filtering: ", percent_of_edges_filter

        tot_num_genes = len(ranking['dest'].unique())
        genes = ranking['dest'].unique()
        tot_num_regulated_edges = 0
        tot_num_validated_edges = 0



        TFs_in_Network = ranking["source"].unique()

        validated_edges_df = pd.DataFrame(columns=["Source Gene","Destination Gene", "Validation_type", "FoldChange/FDR", "Score", "Source Description", "Destination Description"])

        #Validate TARGET
        targetspath="./TFTargets/"

        num_tfs_to_validate = 0

        for file in os.listdir(targetspath):
            if file.startswith("AT"):
                parts = file.split("_")
                tf_name = parts[0]
                if tf_name in TFs_in_Network and tf_name in tfs_list_tovalidate:
                    num_tfs_to_validate = num_tfs_to_validate + 1

                    table = pd.read_table(targetspath+file, delimiter = ",")
                    table.columns = ["target", "foldchange", "fdr"]
                    
                    #Filter pvalue less than 0.05
                    table = table[table['fdr']<pvalue_thres]
                    table.reset_index(inplace=True)
                    table.drop("index", axis=1, inplace=True)

                    table["target"] = map(lambda x: x.upper(), table["target"])
                    targets_list = table["target"]
                    net_targets = self.getTargets_fromNetwork_dataframe(ranking, tf_name)
                    genesect = self.interesect_twogenelists(net_targets, targets_list)
                    num_rows = validated_edges_df.shape[0]
                    
                    genesect_univ = self.interesect_twogenelists(genes, targets_list)
                    tot_num_regulated_edges = tot_num_regulated_edges + len(genesect_univ)
                    tot_num_validated_edges = tot_num_validated_edges + len(genesect)

                    if len(genesect) >0:
                        print "Gene", self.get_genedescr(tf_name)
                        print "Running..."
                        # print "Num of targets in the network for tf ", tf_name, " is ", len(net_targets)
                        #print "genesect", genesect
                        for i, gene in enumerate(genesect):
                            #"source1","dest1", "Validation_type1", "pvalue/fdr1", "score_RF1", "source_descr1", "dest_descr1"
                            score = self.getgene_rankingscore(ranking, tf_name, gene)
                            validated_edges_df.loc[num_rows+i]=[tf_name,gene, "TARGET", self.getgene_pvaluefdr(table, gene), score, self.get_genedescr(tf_name), self.get_genedescr(gene)]
                else:
                    if tf_name not in tfs_list_tovalidate:
                        if tf_name in TFs_in_Network:
                            print "Warning, ", tf_name, " not in tfs_list_tovalidate !!!"

        if len(tfs_list_tovalidate) != num_tfs_to_validate:
                print "Error, len(tfs_list_tovalidate) different than  num_tfs_to_validate !!!", len(tfs_list_tovalidate), num_tfs_to_validate
                sys.exit()


        validated_edges_df.sort_values("Score", ascending=0, inplace=True)
        validated_edges_df.reset_index(inplace=True)
        validated_edges_df.drop("index", axis=1, inplace=True)
        #validated_edges_df.to_html("Validated_Edges_withoutDapseq"+name_output_file+".html")
        #validated_edges_df.to_html("Validated_Edges_withoutDapseq_crf4"+name_output_file+".html")
        #validated_edges_df.to_csv("Validated_Edges_withoutDapseq"+name_output_file+".txt", sep="\t")




        try:
            percent_validated_edges = 100*(float(tot_num_validated_edges) / float(tot_num_regulated_edges))
        except:
            print "Error No edges in the gold standard"
            sys.exit()
        # print "percent_validated_edges of total number of regulated edges ", percent_validated_edges
        
        return validated_edges_df



    def print_PR_curve(self, tfs_list_tovalidate, ranking, prcurve_namefile, validated_edges_df, cheat_meth, network_namefile, thres_prec, all_random_ordering, plot_one_random_ordering):
        
        
        ranking4TFs = ranking.copy()
        
        TFs_in_Network = ranking4TFs["source"].unique()
        
        targets_in_Network = ranking4TFs["dest"].unique()
        
        
        print "Num of TFs in the network: ", len(TFs_in_Network)
        print "Num of targets in the network: ", len(targets_in_Network)
        
        
        ranking4TFs.set_index(['source'], inplace=True)
        
        ranking4TFs_= pd.DataFrame()
        for tf in tfs_list_tovalidate:
            if tf in TFs_in_Network:
                try:
                    ranking4TFs.loc[tf].shape[1]
                    ranking4TFs_ = pd.concat([ranking4TFs_, ranking4TFs.loc[tf]])
                except:
                    ranking4TFs_.loc[tf] = ranking4TFs.loc[tf]
         
        
        ranking4TFs_.reset_index(inplace=True)    
        ranking4TFs_.sort_values('score', ascending=False, inplace=True)
        ranking4TFs_.reset_index(inplace=True)
        ranking4TFs_.drop('index', axis=1, inplace=True)
        
        
        validated_edges_df.drop(['Validation_type', 'FoldChange/FDR', 'Source Description', 'Destination Description', 'Score'], axis=1, inplace=True)
        
        validated_edges_df_keys = validated_edges_df.set_index(['Source Gene', 'Destination Gene'])
        validated_edges_df_keys['score'] = 1

        # print ranking4TFs_.shape #scores for gold standard
        # print validated_edges_df.shape #gold_standard

        tfs_gs = validated_edges_df_keys.index.get_level_values(0).unique()
        genes_gs = validated_edges_df_keys.index.get_level_values(1).unique()

        ranking4TFs_filtered0rows0cols = pd.DataFrame()
        if cheat_meth == "yes":
            ranking4TFs_filtered0rows0cols = ranking4TFs_[ranking4TFs_['dest'].isin(genes_gs)]
        elif cheat_meth == "no":
            ranking4TFs_filtered0rows0cols = ranking4TFs_.copy()#Cheating method no!(uncomment this line;)

        # print "Number of targets for TFs to validate ", len(ranking4TFs_['dest'].unique())
        # print "Number of TFs to validate ", len(ranking4TFs_['source'].unique())
        # print "Number of TFs in validated list of edges ", len(validated_edges_df_keys.index.get_level_values(0).unique())
        # print "Number of targets in validated list of edges ", len(validated_edges_df_keys.index.get_level_values(1).unique())

        # print "Number of targets for TFs to validate the subset of genes and TFs that have at least one known interaction", len(ranking4TFs_filtered0rows0cols['dest'].unique())
        # print "Number of TFs to validate over the subset of genes and TFs that have at least one known interaction", len(ranking4TFs_filtered0rows0cols['source'].unique())

        # print "Dimension of dataframe with edges to validate", ranking4TFs_.shape #scores for gold standard
        # print "Dimension of dataframe with VALIDATED edges", validated_edges_df.shape #gold_standard
        # print "Dimension of dataframe with edges to validate the subset of genes and TFs that have at least one known interaction", ranking4TFs_filtered0rows0cols.shape #scores for gold standard


        edges_score = list(ranking4TFs_filtered0rows0cols['score'])

        average_precision = ""

        edges_gs = []

        for ind in range(0, ranking4TFs_filtered0rows0cols.shape[0]):

            sour = ranking4TFs_filtered0rows0cols.iloc[ind]['source']
            dest = ranking4TFs_filtered0rows0cols.iloc[ind]['dest']
            try:
                tmp = float(validated_edges_df_keys.loc[sour,dest])
                edges_gs.append(1)
            except:
                edges_gs.append(0)


        precision, recall, thresholds = precision_recall_curve(edges_gs, edges_score)
        
        average_precision = auc(recall, precision)
        
        # print "Dim precision vector, ", len(precision)

        # print "Dim recall vector, ", len(recall)

        # print "Dim thresholds vector, ", len(thresholds)

        # print "---> Dim edges_gs, ", len(edges_gs)

        # print "---> Dim edges_score", len(edges_score)

        num_valid_edges = validated_edges_df.shape[0]

        print "---> Total number of validated edges: ", num_valid_edges

        rnd_guessing = float(num_valid_edges)/float(len(edges_score))

        print "---> Fraction of total true edges over total number of possible edges (i.e. Random guessing) -- Total number of validated edges(i.e. true positives)  num_valid_edges (", num_valid_edges, ") / len(edges_score)(i.e. total num of possible edges) (", len(edges_score), ") = ", rnd_guessing

        
        #Find the first index for which the precision is bigger or equal than thres_prec
        indx = 0
        for i, p in enumerate(precision):
            if p >= thres_prec and p!=1:
                indx = i + 1
                break
        
        if indx != 0:
         
            indx_plot = indx-1

            edges_score_prec_thres = thresholds[indx-1]

            indx_edges_score_prec_thres = np.sum(np.array(edges_score)>=edges_score_prec_thres)

            # print "debug: indx_edges_score_prec_thres, ", indx_edges_score_prec_thres

            tot_num_tp_and_fp = np.sum(edges_gs[:indx_edges_score_prec_thres])

            # print "debug: tot_num_tp_and_fp, ", tot_num_tp_and_fp

            print " ---> Recall at precision ",thres_prec, ": ", recall[indx-1], 

            print "\n ---> More details: At precision ",thres_prec, " , we have out of ", indx_edges_score_prec_thres, " top-edges (TP+FP): \n recall[indx-1] (", recall[indx-1], ") * validated_edges_df.shape[0] (", num_valid_edges, ") = ", recall[indx-1]*num_valid_edges, "edges True Positives"
            
            print "---> precision, recall, threshold(bigger or equal): ", precision[indx-1], recall[indx-1], thresholds[indx-1]


            if len(edges_score) > len(precision):
                indx = indx + (len(edges_score)-len(precision))

            # print "debug: edges_score[-indx-1] ", edges_score[-indx-1]
            # print "debug: edges_score_prec_thres", edges_score_prec_thres
            
            ranking_filtered_precision_thres = ranking[ranking["score"]>=edges_score_prec_thres]

            num_edges_prec_thres = len(ranking_filtered_precision_thres)

            num_targets_prec_thres = len(ranking_filtered_precision_thres["dest"].unique())

            num_tfs_prec_thres = len(ranking_filtered_precision_thres["source"].unique())

            subname_net = "EdgesScoreThres"+str(round(edges_score_prec_thres,4))+"_"+"Num_edges"+str(num_edges_prec_thres)+"_"+"NumTfs"+str(num_tfs_prec_thres)

            subname_net = subname_net+"_"+"NumTargets"+str(num_targets_prec_thres)

            ranking_filtered_precision_thres.to_csv(network_namefile+"_PrecisThres"+str(thres_prec)+"_"+subname_net+".txt", sep="\t", index=False)
            
            

        #Second (looser) threshold
        thres_prec2 = 0.4#rnd_guessing + 1.4*rnd_guessing#0.40*rnd_guessing

        indx2 = 0
        for i, p in enumerate(precision):
            if p >= thres_prec2 and p!=1:
                indx2 = i + 1
                break
        
        indx2_plot = indx2-1

        if indx2 != 0:
         
            edges_score_prec_thres2 = thresholds[indx2-1]

            indx_edges_score_prec_thres2 = np.sum(np.array(edges_score)>=edges_score_prec_thres2)

            # print "debug: indx_edges_score_prec_thres2, ", indx_edges_score_prec_thres2

            tot_num_tp_and_fp2 = np.sum(edges_gs[:indx_edges_score_prec_thres2])

            # print "debug: tot_num_tp_and_fp2, ", tot_num_tp_and_fp2

            print " ---> Recall2 at precision2 ",thres_prec2, ": ", recall[indx2-1], 

            print "\n ---> More details: At precision2 ",thres_prec2, " , we have out of ", indx_edges_score_prec_thres2, " top-edges (TP+FP): \n recall[indx2-1] (", recall[indx2-1], ") * validated_edges_df.shape[0] (", num_valid_edges, ") = ", recall[indx2-1]*num_valid_edges, "edges True Positives"
            
            print "---> precision2, recall2, threshold2(bigger or equal): ", precision[indx2-1], recall[indx2-1], thresholds[indx2-1]

            if len(edges_score) > len(precision):
                indx2 = indx2 + (len(edges_score)-len(precision))

            # print "debug: edges_score[-indx-1] ", edges_score[-indx2-1]
            # print "debug: edges_score_prec_thres", edges_score_prec_thres2
            
            ranking_filtered_precision_thres2 = ranking[ranking["score"]>=edges_score_prec_thres2]

            num_edges_prec_thres2 = len(ranking_filtered_precision_thres2)

            num_targets_prec_thres2 = len(ranking_filtered_precision_thres2["dest"].unique())

            num_tfs_prec_thres2 = len(ranking_filtered_precision_thres2["source"].unique())

            subname_net2 = "EdgesScoreThres"+str(round(edges_score_prec_thres2,4))+"_"+"Num_edges"+str(num_edges_prec_thres2)+"_"+"NumTfs"+str(num_tfs_prec_thres2)

            subname_net2 = subname_net2+"_"+"NumTargets"+str(num_targets_prec_thres2)

            ranking_filtered_precision_thres2.to_csv(network_namefile+"_PrecisThres"+str(thres_prec2)+"_"+subname_net2+".txt", sep="\t", index=False)
      

        # Plot Precision-Recall curve
        plt.clf()
        plt.close()
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.plot(recall, precision, label='Precision-Recall curve')
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.ylim([0.0, 1.05])
        plt.xlim([0.0, 1.0])
        plt.title('Precision-Recall: AUC={0:0.3f}'.format(average_precision))
        plt.legend(loc="upper right")
        #plt.show()
        plt.savefig(prcurve_namefile+".pdf", format='pdf')

                
        
        #All RANDOM ORDERING PR CURVE
        if all_random_ordering:
            df_dfg = pd.DataFrame(columns=range(2))

            df_dfg[0] = recall
            df_dfg[1] = precision

            df_dfg.to_csv('./results/precision_recall_dfg.txt', sep="\t")    
        
            edges_gs_random = list(edges_gs)

            df_recall = pd.DataFrame(columns=range(1000))
            df_precision = pd.DataFrame(columns=range(1000))

            average_aupr_rnd = []

            rounds = 1000

            count_ = 0
            for i in range(rounds):
                shuffle(edges_gs_random)
                precision_rnd, recall_rnd, _ = precision_recall_curve(edges_gs_random, edges_score)
                average_precision_rnd = auc(recall_rnd, precision_rnd)

                average_aupr_rnd.append(average_precision_rnd)

                df_recall[i] = recall_rnd

                df_precision[i] = precision_rnd

                if average_precision_rnd > average_precision:
                    count_ += 1

            pvalue = float(count_)/float(rounds) 

            print "---> The average AUPR value for random ordering with ", rounds, "rounds ", "is ", np.mean(average_aupr_rnd)

            print "The median AUPR value for random ordering with ", rounds, "rounds ", "is ", np.median(average_aupr_rnd)

            df_recall.to_csv('./results/random_ordering_recalls.txt', sep="\t")
            df_precision.to_csv('./results/random_ordering_precisions.txt', sep="\t")

            print "The pvalue is: ", pvalue

        if plot_one_random_ordering:

            edges_gs_random = list(edges_gs)
            shuffle(edges_gs_random)
            precision_rnd, recall_rnd, _ = precision_recall_curve(edges_gs_random, edges_score)
            average_precision_rnd = auc(recall_rnd, precision_rnd)
            
            # Plot Precision-Recall curve
            plt.clf()
            plt.close()
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            ax1.plot(recall_rnd, precision_rnd, label='Precision-Recall curve')
            plt.xlabel('Recall')
            plt.ylabel('Precision')
            plt.ylim([0.0, 1.05])
            plt.xlim([0.0, 1.0])
            if all_random_ordering:
                plt.title('Precision-Recall Rnd Ordering: AUC={0:0.3f}'.format(average_precision_rnd) + ' - Pvalue={0:0.2f}'.format(pvalue))
            else:
                plt.title('Precision-Recall Rnd Ordering: AUC={0:0.3f}'.format(average_precision_rnd))

            plt.legend(loc="upper right")
            prcurve_namefile = prcurve_namefile+"RandomOrdering"
            plt.savefig(prcurve_namefile+".pdf", format='pdf')


            # Plot Precision-Recall curve
            plt.clf()
            plt.close()
            fig = plt.figure()
            ax1 = fig.add_subplot(111)

            ax1.plot(recall, precision)

            if indx!=0:
                ax1.plot(recall[indx_plot], precision[indx_plot], 'ro')
                prcurve_namefile = prcurve_namefile+"_Prec1_{0:0.4f}".format(precision[indx_plot])


            ax1.plot(recall_rnd, precision_rnd)

            plt.xlabel('Recall')
            plt.ylabel('Precision')
            plt.ylim([0.0, 1.05])
            plt.xlim([0.0, 1.0])
            plt.title('Precision-Recall Curves')
            plt.legend(loc="upper right")
            prcurve_namefile = prcurve_namefile+"Overlapped"
            plt.savefig(prcurve_namefile+".pdf", format='pdf')
            print "Pipeline End"
        
        return validated_edges_df_keys, ranking4TFs_filtered0rows0cols, edges_gs, edges_score, precision, recall, average_precision







