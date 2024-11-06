import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
# from exp.flux_dep_qubit_spec_J import *
# raw_data = np.load(r'q2_t1_raw_20240129-201538.npz', allow_pickle=True)# ["arr_0"].item()
# tomo_data =
# for k, v in raw_data.items():
#     print(k, v.shape)
# threshold = [2.007e-04, -5.748e-06, 1.421e-04]
#a = raw_data["q2_ro"]
# qf=[4.1985, 4.1965, 4.1965, 4.1985, 4.1985, 4.194500000000001, 4.1925, 4.1925, 4.1865000000000006, 4.1845, 
#     4.1805, 4.1765, 3.9385000000000003, 4.1685, 4.1625000000000005, 4.1585, 4.1525, 4.1465000000000005, 
#     4.1405, 4.1325, 4.1265, 4.1185, 4.1105, 4.1025, 4.0945, 4.0845, 4.0745000000000005, 4.0725, 4.0645, 
#     4.0485, 4.0365, 4.0265, 4.0145, 4.0005, 3.9885, 3.9785, 3.9645, 3.9525, 3.9385000000000003, 3.9245, 3.9105]
#for i in range(50,91):
    #index = np.argmin(a[1,i,0:151])
   # qf.append(4.1985-0.002*(150-index))
#print(qf)
#plt.plot(dfs,a[1,53,:])

from qualang_tools.plot.fitting import Fit
def fit_T1( evo_time, signal ):
    fit = Fit()
    decay_fit = fit.T1( evo_time, signal )
    relaxation_time = np.round(np.abs(decay_fit["T1"][0]) / 4) * 4
    fit_func = decay_fit["fit_func"]
    return relaxation_time, fit_func

def plot_T1( x, y, y_label:list=["I","Q"], fig=None ):
    """

    x shape (M,) 1D array
    y shape (N,M)
    N is 1(I only) or 2(both IQ)
    """
    signal_num = 2
    if fig == None:
        fig, ax = plt.subplots(nrows=signal_num)
    # c = ax.pcolormesh(dfs, amp_log_ratio, np.abs(s21), cmap='RdBu')# , vmin=z_min, vmax=z_max)
    # ax.set_title('pcolormesh')
    # fig.show()
    # Plot
    fig.suptitle("T1 measurement")
    for i in range(signal_num):
        ax[i].plot( x, y, label="data")
        ax[i].set_ylabel(f"{y_label} quadrature [V]")
        ax[i].set_xlabel("Wait time (ns)")

        fit_T1_par, fit_func = fit_T1(x, y)
        ax[i].plot( x, fit_func(x), label="fit")

    return fig

def anal_t1(raw_data, ro_name, t_delay, repeat_time):
    """
    raw data shape (repeat time, IQ, t_delay)

    """
    t1_collection = []
    for i in range(repeat_time):
        iqvector = raw_data[ro_name][i][0][0]+1j*raw_data[ro_name][i][1][0]-(raw_data[ro_name][i][0][-1]+1j*raw_data[ro_name][i][1][-1])
        ivector = raw_data[ro_name][i][0][0] - raw_data[ro_name][i][0][-1]
        theta = np.arccos(np.abs(ivector)/np.abs(iqvector))
        
        new_iqdata = raw_data[ro_name][i][0]+1j*raw_data[ro_name][i][1]-(raw_data[ro_name][i][0][-1]+1j*raw_data[ro_name][i][1][-1])
        rotate_iqdata = np.exp(1j*theta)*new_iqdata
        new_iqdata_2 = raw_data[ro_name][i][0][-1]+1j*raw_data[ro_name][i][1][-1] + rotate_iqdata
        iqdata_ex = raw_data[ro_name][i][0] +1j*raw_data[ro_name][i][1]
        #print(new_iqdata)
        ampdata = np.real(new_iqdata_2)
        T1_i = fit_T1(t_delay*4, ampdata)[0]
        t1_collection.append(T1_i)
        #plot_T1(t_delay,ampdata)
    t1_collection = np.array(t1_collection)/1000
    return t1_collection

def anal_tune_flux_t1(raw_data, ro_name, t_delay, flux,repeat_time):
    """

    raw data shape (flux, repeat time, IQ, t_delay)

    """
    t1_collection = []
    t1_ave=[]
    for i in range(len(flux)):
        t1_collection.append([])
    for j in range(len(flux)):
        a=[]
        for i in range(repeat_time):
            iqvector = raw_data[ro_name][j][i][0][0]+1j*raw_data[ro_name][j][i][1][0]-(raw_data[ro_name][j][i][0][-1]+1j*raw_data[ro_name][j][i][1][-1])
            ivector = raw_data[ro_name][j][i][0][0] - raw_data[ro_name][j][i][0][-1]
            theta = np.arccos(np.abs(ivector)/np.abs(iqvector))
            
            new_iqdata = raw_data[ro_name][j][i][0]+1j*raw_data[ro_name][j][i][1]-(raw_data[ro_name][j][i][0][-1]+1j*raw_data[ro_name][j][i][1][-1])
            rotate_iqdata = np.exp(1j*theta)*new_iqdata
            new_iqdata_2 = raw_data[ro_name][j][i][0][-1]+1j*raw_data[ro_name][j][i][1][-1] + rotate_iqdata
            iqdata_ex = raw_data[ro_name][j][i][0] +1j*raw_data[ro_name][j][i][1]
            #print(new_iqdata)
            ampdata = np.real(new_iqdata_2)
            a.append(ampdata)
            T1_i = fit_T1(t_delay*4, ampdata)[0]
            t1_collection[j].append(T1_i)
        a = np.array(a)
        b = np.mean(a, axis=0)
        t1_ave.append(b)
    t1_ave = np.array(t1_ave)
    print(t1_ave.shape)
    t1_collection = np.array(t1_collection)/1000
    t1_collection = tune_flux_t1_correction(t1_collection,flux,repeat_time)
    return t1_collection, t1_ave

def gaussian_pdf_cdf(x,mu,sigma):
    pdf = 1/(sigma * np.sqrt(2 * np.pi))*np.exp( - (x - mu)**2 / (2 * sigma**2) )
    cdf = pdf.cumsum()
    cdf /= cdf[-1]
    return pdf, cdf
def t1_hist(data):
    fig, ax = plt.subplots()
    #bin_width = 100
    #start_value = 0
    #end_value = 900
    #custom_bins = [start_value + i * bin_width for i in range(int((end_value - start_value) / bin_width) + 1)]
    #hist_values, bin_edges = np.histogram(a, bins=custom_bins, density=True)
    data_norm = (data-np.mean(data))/np.std(data)+np.mean(data)
    t1_mean = np.mean(data)
    sigma = np.std(data)
    print(np.mean(data))
    count, bins, ignore=ax.hist(data, 30,density=True, cumulative=True,alpha=0.8, color='blue', label='Histogram_CDF')
    ax.set_xlabel("T1 (us)")
    ax.set_ylabel("cdf")
    for i in range(len(count)):
        if count[i] >= 0.5:
            print("median=%.2f"%bins[i])
            break
    xmin, xmax = ax.get_xlim()
    pdf,cdf = gaussian_pdf_cdf(bins,t1_mean,sigma)
    plt.plot(bins, pdf, linewidth=2, color='r',label='PDF_fit')
    plt.plot(bins, cdf, color='y', label='CDF_fit')
    fig.suptitle('T1 Distribution')

def tune_flux_t1_correction(data,flux,repeat_time):
    t1=[]
    for i in range(len(flux)):
        b=[]
        for j in range(repeat_time):
            if data[i][j]<=35:
                b.append(data[i][j])
        #b=np.array(b)
        t1.append(b)
    t1 = np.array(t1)
    return t1

def tune_flux_t1_spectrum(data,flux,repeat_time):
    T1_avg=[]
    error=[]
    for i in range(len(flux)):
        T1_avg.append(np.mean(data[i]))
        error.append(np.max(data[i])-np.mean(data[i]))
        print(np.mean(data[i]))
    plt.errorbar(flux,T1_avg,error,elinewidth=2,capsize=4,fmt='o')

def plot_qubit_flux_decay( data, flux, t_delay):
    """
    data shape ( flux, t_delay )

    """
    fig, ax = plt.subplots()
    ax.set_title('pcolormesh')
    ax.set_xlabel("T1 (us)")
    ax.set_ylabel("Flux")
    ax.pcolormesh( t_delay/1000, flux, data, cmap='RdBu')# , vmin=z_min, vmax=z_max)
    T1 = []
    for i in range(len(flux)):
        try:
            T1_i = fit_T1(t_delay, data[i,:])[0]/1000
        except:
            T1_i = 0
        T1.append(T1_i)
    T1 = np.array(T1)
    print(T1)
    ax.plot(T1,flux, color='b',linewidth=1.5, label="T1")
    ax.set_xlim(0,0.4)
    #ax.set_xlim([0,20])

if __name__ == '__main__':
    dataset = xr.open_dataset("20240625_1425_q0_xy_T1spectrum.nc")
    for ro_name, data in dataset.data_vars.items():
        idata = data.values[0]
    #pan = 300 * u.MHz
    #df = 2 * u.MHz
    flux_range = (-0.3, 0.3)
    flux_resolu = 0.003
    flux_idle = -0.3473
    tau_min = 16 // 4 # in clock cycles
    tau_max = 80000 // 4  # in clock cycles
    d_tau = 1600 // 4  # in clock cycles

    ro_name = "q2_ro"
    repeat_time = 10

    t_delay = np.arange(tau_min, tau_max + 0.1, d_tau)  # Linear sweep
    time = dataset.coords["time"].values
    flux = np.arange(flux_range[0],flux_range[1],flux_resolu)
    #dfs = np.arange(-span, +span, df)

    ########################## t1 data ################################################
    #data = anal_t1(raw_data,ro_name,t_delay,repeat_time)

    ########################## flux-t1 data ###########################################
    #data, anal_raw_data = anal_tune_flux_t1(raw_data,ro_name,t_delay, flux, repeat_time)
    
    ########################## t1 cdf histogram #######################################
    #t1_hist(data)

    ########################## flux-t1 histogram ######################################
    #for i in range(len(flux)):
        #print(i)
        #t1_hist(t1_data[i])
    
    ########################## flux-t1 spectrum #######################################
    #tune_flux_t1_spectrum(data,flux,repeat_time)
    
    ########################## flux-t1 mix cdf histogram ##############################
    #t1_mix=[]
    #for i in range(len(flux)):
        #t1_mix.extend(list(data[i]))
    #t1_mix=np.array(t1_mix)
    #t1_hist(t1_mix)

    ############################### plot flux-t1 ######################################
    #t1_data_ave = []
    #for i in range(len(flux)):
        #t1_data_ave.append(np.mean(data[i]))
    #t1_data_ave = np.array(t1_data_ave)
    plot_qubit_flux_decay(idata,flux,time)

    ##################################################################################
    plt.legend()
    plt.show()
