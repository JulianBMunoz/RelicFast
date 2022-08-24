color=sns.color_palette('flare', 7)

for idx in range(100):

    lcdm = np.loadtxt("/Users/nicholasdeporzio/Downloads/test_lcdm/_matterpower_"+str(idx+1)+".dat")
    plusde1 = np.loadtxt("/Users/nicholasdeporzio/Downloads/test_plusde1/_matterpower_"+str(idx+1)+".dat")
    plusde2 = np.loadtxt("/Users/nicholasdeporzio/Downloads/test_plusde2/_matterpower_"+str(idx+1)+".dat")
    plusde3 = np.loadtxt("/Users/nicholasdeporzio/Downloads/test_plusde3/_matterpower_"+str(idx+1)+".dat")
    plusde4 = np.loadtxt("/Users/nicholasdeporzio/Downloads/test_plusde4/_matterpower_"+str(idx+1)+".dat")
    plusde5 = np.loadtxt("/Users/nicholasdeporzio/Downloads/test_plusde5/_matterpower_"+str(idx+1)+".dat")
    plusde6 = np.loadtxt("/Users/nicholasdeporzio/Downloads/test_plusde6/_matterpower_"+str(idx+1)+".dat")
    
    lcdm_interp = scipy.interpolate.interp1d(lcdm[:,0]*0.70148, lcdm[:,1])
    plusde_interp1 = scipy.interpolate.interp1d(plusde1[:,0]*0.70148, plusde1[:,1])
    plusde_interp2 = scipy.interpolate.interp1d(plusde2[:,0]*0.70148, plusde2[:,1])
    plusde_interp3 = scipy.interpolate.interp1d(plusde3[:,0]*0.70148, plusde3[:,1])
    plusde_interp4 = scipy.interpolate.interp1d(plusde4[:,0]*0.70148, plusde4[:,1])
    plusde_interp5 = scipy.interpolate.interp1d(plusde5[:,0]*0.70148, plusde5[:,1])
    plusde_interp6 = scipy.interpolate.interp1d(plusde6[:,0]*0.70148, plusde6[:,1])
    
    plt.figure(figsize=(15, 7.5))
    plt.title("CAMB Output, z="+f"{z_vals[100-idx-1]:.3f}"+", Normalized at Small k")
    plt.plot(plot_x, lcdm_interp(plot_x)/lcdm_interp(plot_x[0]), label='LCDM', color=colors[0])
    plt.plot(plot_x, plusde_interp1(plot_x)/plusde_interp1(plot_x[0]), label="5%+DE", color=colors[1])
    plt.plot(plot_x, plusde_interp2(plot_x)/plusde_interp2(plot_x[0]), label="10%+DE", color=colors[2])
    plt.plot(plot_x, plusde_interp3(plot_x)/plusde_interp3(plot_x[0]), label="15%+DE", color=colors[3])
    plt.plot(plot_x, plusde_interp4(plot_x)/plusde_interp4(plot_x[0]), label="20%+DE", color=colors[4])
    plt.plot(plot_x, plusde_interp5(plot_x)/plusde_interp5(plot_x[0]), label="25%+DE", color=colors[5])
    plt.plot(plot_x, plusde_interp6(plot_x)/plusde_interp6(plot_x[0]), label="30%+DE", color=colors[6])
    plt.plot([0.05, 0.05], [min(lcdm_interp(plot_x)/lcdm_interp(plot_x[0])), max(lcdm_interp(plot_x)/lcdm_interp(plot_x[0]))], label="$k_p$", color='black')
    plt.plot([0.2*0.70148, 0.2*0.70148], [min(lcdm_interp(plot_x)/lcdm_interp(plot_x[0])), max(lcdm_interp(plot_x)/lcdm_interp(plot_x[0]))], label="$k_{nl}$", color='black', linestyle='dashed')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('k')
    plt.ylabel('Pmm')
    plt.legend()
    plt.savefig("/Users/nicholasdeporzio/Desktop/plots/pmm/pmm_"+str(idx+1)+".png")