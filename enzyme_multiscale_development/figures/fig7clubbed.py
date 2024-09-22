#!/usr/bin/env python3

from figsurfplus import *

if __name__ == '__main__':
    aminos_3letter = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS',
                        'ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP',
                        'TYR','VAL',]
    aminos_1letter = [ i for i in 'ARNDCQEGHILKMFPSTWYV' ]
    aminos_switch_symbol = dict(zip(aminos_3letter,aminos_1letter)) | \
        dict(zip(aminos_1letter,aminos_3letter))
    
    hydrophobics_kd = set([ i for i in 'IVLFCMA' ])
    positives = set([ i for i in 'KRH' ])
    negatives = set([ i for i in 'DE' ])
#    others = set(aminos_1letter) - hydrophobics_kd - positives - negatives
    others = set([ i for i in 'KRHDE' ]) # testing 0716
    
    proteins = ['lip','dha']
    proteinnames = ['Lipase','Dehalogenase']
    
    resolutions = ['aa','cg']
    resolutionnames = ['AA','CG']
    
    data = [
        'data/rtypes_counts_ensemble-lip-aa.csv',
        'data/rtypes_counts_ensemble-lip-cg.csv',
        'data/rtypes_counts_ensemble-dha-aa.csv',
        'data/rtypes_counts_ensemble-dha-cg.csv',
    ]
    
    
    patchlist = []
    amino_colors = ['orange', 'blue', 'black', 'red', 'black', \
        'black', 'red', 'black', 'blue', 'orange', \
        'orange', 'blue', 'orange', 'orange', 'black', \
        'black', 'black', 'orange', 'orange', 'orange']
    for subsetno,subset in enumerate([hydrophobics_kd, positives, negatives, others]):
        fig7c,ax7c = plt.subplots(1,1, figsize=(4.5,9/16*4.5))

        rePaired = sns.color_palette('Paired')[0:2]+sns.color_palette('Paired')[4:6]
        ccax5 = (cycler(color=rePaired) +
         cycler(ls=['--', '-']*2))
        ax7c.set_prop_cycle(ccax5)
        
        for filepath in data:
            mydata = np.loadtxt(filepath, comments=('@','#'), delimiter=',')
            if 'aa' in filepath: myres='AA'
            if 'cg' in filepath: myres='CG'
            if 'aa' in filepath:
                ls='--'
                res='AA'
            if 'cg' in filepath:
                ls='-'
                res='CG'  
            if 'lip' in filepath: sysname = 'Lipase'
            if 'dha' in filepath: sysname = 'Dehalogenase'

            reduced_list = []
            reduced_data = np.full(mydata[:,0].shape,None)
            
            for idx,amino in enumerate(aminos_1letter):
                if amino in subset:
                    reduced_list.append(aminos_switch_symbol[amino])
                    if reduced_data.any():
                        reduced_data = np.column_stack((reduced_data,mydata[:,idx]))
                    else: reduced_data = mydata[:,idx]
            print(reduced_list)
    #        print(reduced_data)
            
    #        patchname = plot_to_ax(reduced_list,reduced_data,ax7c,'violin',
    #            labels=reduced_list,resolution=myres)
    #        patchlist.append(patchname)

            surfres_count_ensemble = np.sum(mydata,axis=1)
            type_proportion = reduced_data / surfres_count_ensemble[:,None]
            mean_type_proportion = np.mean(type_proportion,axis=0)
            mean_type_stderr = np.std(type_proportion,axis=0)/np.sqrt(len(type_proportion))
        #                axexp2.plot(amino_acid_list,mean_type_proportion,linestyle='None')
            ax7c.errorbar(reduced_list, mean_type_proportion, yerr=2*1.96*mean_type_stderr,linestyle='None',marker='.',label=f'{sysname} {res}')

#        if subsetno==0: ax7c.legend()
        if subsetno in [3]: ax7c.legend(borderpad=0.2, labelspacing=0.2 )
#        ax7c.set_title(f'Mean Proportion of Total Surface Residues')#+\
#                        f'distribution of mean surface counts')
#        ax7c.set_xlabel(f'amino acid')
#        ax7c.set_ylabel(f'surface count')                        
#            
        plt.tight_layout()       
        fig7c.savefig(f'images/fig7.png', dpi=300)

        plt.close()        


