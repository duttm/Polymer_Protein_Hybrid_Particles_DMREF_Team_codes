#!/usr/bin/env python3

from fig2plus import *
mpl.rcParams.update({'font.size': 16})


if __name__ == '__main__':
    params = ['circular_variance','kyte_doolittle','electrostatics']
    paramnames = ['Circular Variance','Hydrophobicity Index','Electrostatic Potential [kT/e]']
    xlims = [[0,1],[-4.5,4.5],[-20,20]]
    
    paramprint = dict(zip(params,paramnames))
    paramrange = dict(zip(params,xlims))
    
    proteins = ['lip','dha']
    proteinnames = ['Lipase','Dehalogenase']
    
    resolutions = ['aa','cg']
    resolutionnames = ['AA','CG']
    
    for param in params:
    
        data = [
            'data/paramvals-combined-'+param+'-lip-aa.csv',
            'data/paramvals-combined-'+param+'-lip-cg.csv',
            'data/paramvals-combined-'+param+'-dha-aa.csv',
            'data/paramvals-combined-'+param+'-dha-cg.csv'
        ]
        
        fig5,ax5 = plt.subplots(1,1, figsize=(4,4))
#        fig5,ax5 = plt.subplots(1,1, figsize=(4.5,9/16*4.5))
        rePaired = sns.color_palette('Paired')[0:2]+sns.color_palette('Paired')[4:6]
        ccax5 = (cycler(color=rePaired) +
         cycler(ls=['--', '-']*2))
        ax5.set_prop_cycle(ccax5)

        for filepath in data:
            mydata = np.loadtxt(filepath, comments=('@','#'))
            bins=50
            if 'elec' in filepath: bins = 100
            y,binEdges = np.histogram(mydata,bins=bins,density=True)
            myerror = 0
            if 'aa' in filepath:
                ls=':'
                res='AA'
            if 'cg' in filepath:
                ls='-'
                res='CG'  
            if 'lip' in filepath: sysname = 'Lip'
            if 'dha' in filepath: sysname = 'Dha'
            if paramrange[param]: ax5.set_xlim(paramrange[param])
            
            ax5.stairs(y,binEdges,label=f'{sysname} {res}',linestyle=ls,lw=3)
            
#            ax5.legend(frameon=False, borderpad=0, handlelength=1.5, handletextpad=0.5)
            if param in ['kyte_doolittle']: ax5.legend(frameon=False, borderpad=0)

            
            
#            ax5.set_title(f'Surface {paramprint[param]} - AA vs CG')#\n'+\
#                            +f'pooled ensemble histogram')
            ax5.set_xlabel(f'{paramprint[param]}')
                
            plt.tight_layout()       
            fig5.savefig(f'images/fig5-{param}.png', dpi=300)

            plt.close()        


