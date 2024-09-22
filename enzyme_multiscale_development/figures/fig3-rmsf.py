import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import t
mpl.rcParams.update({'font.size': 15})


def read_file(file_path):
    # Function to read a file into a numpy array with 2 columns
    data = np.loadtxt(file_path, comments=("@","#"))
    return data

def combine_data(file_paths):
    # Function to combine the second column from each array into a single dataset
    combined_data = np.column_stack([read_file(file_path)[:, 1] for file_path in file_paths])
    return combined_data

def compute_confidence_interval(data, confidence_level=0.95):
    # Function to compute confidence interval for each row
    n = data.shape[1]
    means = np.mean(data, axis=1)
    std_devs = np.std(data, axis=1, ddof=1)
    margin_of_errors = t.ppf((1 + confidence_level) / 2, n - 1) * (std_devs / np.sqrt(n))
    confidence_intervals = np.column_stack((means - margin_of_errors, means + margin_of_errors))
    return means,confidence_intervals

def plot_confidence_intervals(means,confidence_intervals):
    # Function to plot confidence intervals
    x_values = np.arange(1, len(confidence_intervals) + 1)
#    plt.plot(x_values, means, label='Mean', linestyle="-")
    plt.plot(x_values, means, label='Mean',lw=2.5)
#    plt.fill_between(x_values, confidence_intervals[:, 0], confidence_intervals[:, 1], alpha=0.3, label='95% CI')

#    plt.xlabel('residue number')
#    plt.ylabel('RMSF Mean Values with 95% CI')
#    plt.title('95% Confidence Intervals for Mean Values by Row')
#    plt.legend()
#    plt.savefig('conf.png')
#    plt.show()

# Example usage
plt.cla()

#file_paths = ['/home/mason/dmref/lipase/aa/production/traj_samples/analysis_data/rmsf/rmsf-seed'+str(i)+'-Calpha.xvg' for i in range (1,11)]  # Replace with your actual file paths
#aacombined_data = combine_data(file_paths)
#means,confidence_intervals = compute_confidence_interval(aacombined_data)

#plot_confidence_intervals(means,confidence_intervals)

#file_paths = ['/home/mason/dmref/lipase/cg/production/traj_samples/analysis_data/rmsf/rmsf-seed'+str(i)+'-BB.xvg' for i in range (1,11)]  # Replace with your actual file paths
#cgcombined_data = combine_data(file_paths)
#means,confidence_intervals = compute_confidence_interval(cgcombined_data)

#plot_confidence_intervals(means,confidence_intervals)


my_sys = 'dehalogenase'

if my_sys=='lipase':
    file_paths = ['/home/mason/dmref/lipase/aa/production/analysis_data/rmsf/rmsf-seed'+str(i)+'-Calpha.xvg' for i in range (1,11)]  # Replace with your actual file paths
if my_sys=='dehalogenase':
    file_paths = ['/home/mason/dmref/dehalogenase/aa/production/analysis_data/rmsf/rmsf-seed'+str(i)+'-Calpha.xvg' for i in range (1,11)]  # Replace with your actual file paths

aacombined_data = combine_data(file_paths)
means,confidence_intervals = compute_confidence_interval(aacombined_data)

plot_confidence_intervals(means,confidence_intervals)

if my_sys=='lipase':
    file_paths = ['/home/mason/dmref/lipase/cg/production/analysis_data/rmsf/rmsf-seed'+str(i)+'-BB.xvg' for i in range (1,11)]  # Replace with your actual file paths
if my_sys=='dehalogenase':
    file_paths = ['/home/mason/dmref/dehalogenase/cg/production/analysis_data/rmsf/rmsf-seed'+str(i)+'-BB.xvg' for i in range (1,11)]  # Replace with your actual file paths

cgcombined_data = combine_data(file_paths)
means,confidence_intervals = compute_confidence_interval(cgcombined_data)

plot_confidence_intervals(means,confidence_intervals)

#file_paths = ['/home/mason/dmref/figs/1oil-pdb-calpha-tempfactors.txt']  # Replace with your actual file paths
#cgcombined_data = combine_data(file_paths)
if my_sys=='lipase':
    tempfactordata = read_file('/home/mason/dmref/figs/1oil-pdb-calpha-tempfactors.txt')
if my_sys=='dehalogenase':
    tempfactordata = read_file('/home/mason/dmref/figs/3rk4-pdb-calpha-tempfactors.txt')
#B = np.reshape(tempfactordata,(-1,2))
#means,confidence_intervals = compute_confidence_interval(B)
plt.plot(tempfactordata/np.max(tempfactordata)/2,linestyle=(0,(1,1)), linewidth=3)
#plot_confidence_intervals(means,confidence_intervals)



plt.xlabel('Residue Number')
plt.ylabel('RMSF')
#plt.legend(f'{my_sys}')
#plt.legend(['AA','CG','PDB'],loc='upper center', borderpad=0.3, labelspacing=0.5, frameon=True, borderaxespad=-3.5, handletextpad=0.5,ncol=3, fancybox=True, shadow=True)
fig = plt.gcf()
fig.set_size_inches(4.5, 3/4*4.5)
plt.tight_layout()
plt.savefig('images/fig5-'+my_sys+'-testing0720.png',dpi=300)
#plt.show()

