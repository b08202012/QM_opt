import numpy as np
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture


class GMM_model():

    def __init__( self ):
        self.gmm = GaussianMixture(n_components=2, random_state=0)

    def import_trainingData( self, data ):
        """
        input numpy array with shape (n,2)
        n is point number
        """
        self.training_data = data
        self.gmm.fit(data)
        self.gmm.weights_ = [0.5,0.5]

        # return self
    def relabel_model( self, ground_data ):
        """
        input numpy array with shape (2,n)
        """
        print( ground_data.shape )

        gp = np.array([np.mean(ground_data, axis=1)])
        # print( gp )
        # print( gp.shape )

        # print(self.gmm.predict( gp ))
        if self.gmm.predict( gp ) == 1:
            self.gmm.means_ = np.flip(self.gmm.means_,0)
            self.gmm.weights_ = np.flip(self.gmm.weights_,0)
            self.gmm.covariances_ = np.flip(self.gmm.covariances_,0)
            self.gmm.precisions_cholesky_ = np.flip(self.gmm.precisions_cholesky_,0)

    def output_paras( self ):
        """
        four para in dict
        means
        weights
        covariances
        precisions_cholesk
        """
        output_dict = {
            "means":self.gmm.means_,
            "weights":self.gmm.weights_,
            "covariances":self.gmm.covariances_,
            "precisions_cholesky":self.gmm.precisions_cholesky_,
        }

        return output_dict

    def rebuild_model( self, paras:dict ):

        self.gmm.means_ = paras["means"]
        self.gmm.weights_ = paras["weights"]
        self.gmm.covariances_ = paras["covariances"]
        self.gmm.precisions_cholesky_ = paras["precisions_cholesky"]


    def get_prediction( self, data ):
        """
        input numpy array with shape (n,2)
        n is point number
        """
        self.__input_data = data
        self.__predict_label = self.gmm.predict( data )
        return self.__predict_label

    def get_label( self ):
        """
        input numpy array with shape (n,2)
        n is point number
        """

        return self.__predict_label

    def get_distribution( self ):
        return np.bincount(self.__predict_label)

def train_model( data ):
    """
    data type
    3 dim with shape (2*2*N)
    shape[0] is prepare state
    shape[1] is I and Q
    shape[2] is N times single shot
    
    """
    idata = np.array([])
    qdata = np.array([])
    for i in range(data.shape[0]):
        idata = np.append(idata,data[i][0])
        qdata = np.append(qdata,data[i][1])
    training_data = np.array([idata, qdata])
    my_model = GMM_model()
    my_model.import_trainingData(training_data.transpose())
    my_model.relabel_model(data[0])
    return my_model

def create_img( data, dist_model, fig=None, label="" ):

    """
    data type
    3 dim with shape (2*2*N)
    shape[0] is prepare state
    shape[1] is I and Q
    shape[2] is N times single shot
    
    """

    plt.rcParams['figure.figsize'] = [16.0, 8.0]
    plt.rcParams['figure.autolayout'] = True
    if fig==None:
        fig = plt.figure(constrained_layout=True)
    gs = fig.add_gridspec(2, 3)
    ax_iq_all = plt.subplot(gs[0, 0])
    ax_iq_all.set_title(f"{label} Training data", fontsize=20  )
    ax_iq_all.tick_params(axis='both', labelsize=15)
    # ax_iq_all.tick_params(axis='y', labelsize=15)

    # ax_iq_all.yticks(fontsize=20)
    ax_iq_0 = plt.subplot(gs[0,1])
    ax_iq_0.set_title("Prepare |0>", fontsize=20  )
    ax_iq_0.tick_params(axis='both', labelsize=15)
    ax_hist_0 = plt.subplot(gs[1,1])
    ax_hist_0.tick_params(axis='both', labelsize=15)

    ax_iq_1 = plt.subplot(gs[0,2])
    ax_iq_1.set_title("Prepare |1>", fontsize=20  )
    ax_iq_1.tick_params(axis='both', labelsize=15)

    ax_hist_1 = plt.subplot(gs[1,2])
    ax_hist_1.tick_params(axis='both', labelsize=15)

    # scatter plot
    training_data = dist_model.training_data.transpose()
    training_plot = {
        "I":training_data[0],
        "Q":training_data[1],
        "label": dist_model.get_prediction(dist_model.training_data)
    }
    make_scatter_dist( training_plot, ax_iq_all )

    make_ellipses( dist_model.gmm, ax_iq_all )

    prepare_0_data = data[0]
    prepare_0_plot = {
        "I":prepare_0_data[0],
        "Q":prepare_0_data[1],
        "label": dist_model.get_prediction(prepare_0_data.transpose())
    }
    make_scatter_dist(prepare_0_plot,ax_iq_0)
    prepare_0_dist = dist_model.get_distribution()
    prepare_0_dist = prepare_0_dist/np.sum(prepare_0_dist)

    ax_iq_0.text(0.07,0.9,f"P(0|0)={prepare_0_dist[0]:.3f}", fontsize = 20, transform=ax_iq_0.transAxes)
    ax_iq_0.text(0.07,0.8,f"P(0|1)={prepare_0_dist[1]:.3f}", fontsize = 20, transform=ax_iq_0.transAxes)

    prepare_1_data = data[1]
    prepare_1_plot = {
        "I":prepare_1_data[0],
        "Q":prepare_1_data[1],
        "label": dist_model.get_prediction(prepare_1_data.transpose())
    }
    make_scatter_dist(prepare_1_plot,ax_iq_1)
    prepare_1_dist = dist_model.get_distribution()
    prepare_1_dist = prepare_1_dist/np.sum(prepare_1_dist)
    ax_iq_1.text(0.07,0.9,f"P(1|0)={prepare_1_dist[0]:.3f}", fontsize = 20, transform=ax_iq_1.transAxes)
    ax_iq_1.text(0.07,0.8,f"P(1|1)={prepare_1_dist[1]:.3f}", fontsize = 20, transform=ax_iq_1.transAxes)

    # Histogram plot
    sigma_0 = get_sigma(dist_model.output_paras()["covariances"][0])
    sigma_1 = get_sigma(dist_model.output_paras()["covariances"][1])
    sigma = np.max( np.array([sigma_0,sigma_1]) )
    centers = dist_model.output_paras()["means"]
    print(centers)
    pos = get_proj_distance(centers,centers.transpose())
    dis = np.abs(pos[1]-pos[0])
    make_distribution( pos, sigma, get_proj_distance(centers,prepare_0_data), 0, ax_hist_0)
    make_distribution( pos, sigma, get_proj_distance(centers,prepare_1_data), 1, ax_hist_1)

    snr = dis/sigma
    fig.text(0.05,0.35,f"Readout Fidelity={1-(prepare_0_dist[1]+prepare_1_dist[0])/2:.3f}", fontsize = 20)
    fig.text(0.05,0.3,f"IQ distance/STD={dis:.2f}/{sigma:.2f}", fontsize = 20)
    fig.text(0.05,0.25,f"Voltage SNR={snr:.2f}", fontsize = 20)
    fig.text(0.05,0.20,f"Power SNR={np.log10(snr)*20:.2f} dB", fontsize = 20)

    # # Output image
    # img_buf = io.BytesIO()
    # fig.savefig(img_buf, format='png')
    # plt.close(fig)
    return fig

import matplotlib as mpl
colors = ["blue", "red"]

def make_ellipses(gmm, ax):
    for n, color in enumerate(colors):
        match gmm.covariance_type:
            case "full":
                covariances = gmm.covariances_[n][:2, :2]
            case "tied":
                covariances = gmm.covariances_[:2, :2]
            case "diag":
                covariances = np.diag(gmm.covariances_[n][:2])
            case "spherical":
                covariances = np.eye(gmm.means_.shape[1]) * gmm.covariances_[n]
            case _:
                covariances = gmm.covariances_[n][:2, :2]
        v, w = np.linalg.eigh(covariances)
        u = w[0] / np.linalg.norm(w[0])
        angle = np.arctan2(u[1], u[0])
        angle = 180 * angle / np.pi  # convert to degrees
        v = 2.0 * np.sqrt(2.0) * np.sqrt(v)
        ell = mpl.patches.Ellipse(
            gmm.means_[n, :2], v[0], v[1], angle=180 + angle, color=color, fill=False
        )
        ell.set_clip_box(ax.bbox)
        ell.set_alpha(0.5)
        ax.add_artist(ell)
        ax.set_aspect("equal", "datalim")

def get_sigma( covariances ):
    v, w = np.linalg.eigh(covariances)
    v = np.sqrt(v/2)
    return  np.sqrt(v[0]**2+v[1]**2)

def make_scatter_dist(data, ax):
    """
    data type
    3 dim with shape (2*N)
    shape[0] is I and Q
    shape[1] is N times single shot
    
    """
    cmap = mcolors.ListedColormap(["blue", "red"])
    data_pts = data["I"].shape[-1]
    ax.scatter( data["I"], data["Q"], marker='o', c=data["label"], cmap=cmap, s=1000/data_pts, alpha=0.5)
    ax.set_aspect("equal", "datalim")
    ax.set_xlabel('I Voltage Signal', fontsize=20)
    ax.set_ylabel('Q Voltage Signal', fontsize=20)

from lmfit.models import GaussianModel
from lmfit.model import ModelResult
def make_distribution( mu, sigma, data, prepare_state:int, ax):

    dis = np.abs(mu[1]-mu[0])
    est_peak_h = 1/sigma
    # fitting_guess = [1/sigma,pos[0],sigma,1/sigma,pos[1],sigma]
    bin_center = np.linspace(-(dis+2.5*sigma), dis+2.5*sigma,50)

    width = bin_center[1] -bin_center[0]
    bins = np.append(bin_center,bin_center[-1]+width) -width/2

    hist, bin_edges = np.histogram(data, bins, density=True)


    make_histogram(bin_center, hist, ax)
    gm0 = GaussianModel(prefix="g0_", name="g0")
    gm1 = GaussianModel(prefix="g1_", name="g1")
    model = gm0 + gm1
    model.set_param_hint('g0_center',value=mu[0], vary=False)
    model.set_param_hint('g1_center',value=mu[1], vary=False)
    model.set_param_hint('g0_amplitude',min=0, max=est_peak_h*2, vary=True)
    model.set_param_hint('g1_amplitude',min=0, max=est_peak_h*2, vary=True)
    model.set_param_hint('g0_sigma',value=sigma, vary=True)
    model.set_param_hint('g1_sigma',value=sigma, vary=True)

    params = model.make_params()
    results = model.fit(hist,params,x=bin_center)

    # print(results.fit_report())
    
    make_gaussian_fit_curve( bin_center, results, ax)
    ax.set_yscale('log')
    ax.set_ylim(1e-3,est_peak_h*2)
    peak_value = np.array([results.params["g0_amplitude"].value, results.params["g1_amplitude"].value])
    probability = peak_value/np.sum(peak_value)
    ax.text(0.07,0.9,f"P({prepare_state}|0)={probability[0]:.3f}", fontsize = 20, transform=ax.transAxes)
    ax.text(0.07,0.8,f"P({prepare_state}|1)={probability[1]:.3f}", fontsize = 20, transform=ax.transAxes)
    ax.set_xlabel('Projected Voltage Signal', fontsize=20)

def make_histogram(bin_center, data, ax):
    """
    data type
    3 dim with shape (2*N)
    shape[0] is I and Q
    shape[1] is N times single shot
    
    """
    cmap = mcolors.ListedColormap(["blue", "red"])
    width = bin_center[1] -bin_center[0]

    ax.bar(bin_center, data, width=width)


def single_gaussian( x, *params ):
    (a, mu, sigma) = params
    res =   a * np.exp( - (x - mu)**2.0 / (2.0 * sigma**2.0) )
    return res



def make_gaussian_fit_curve(xdata, result:ModelResult, ax):
    # ax.plot(xdata, result.init_fit, '-', label='init fit')     
    ax.plot(xdata, result.best_fit, '--', color="black", label='best fit', linewidth=2)

    # popt, pcov = curve_fit(double_gaussian, xdata, ydata, p0=guess)
    # (c1, mu1, sigma1, c2, mu2, sigma2) = popt
    gm = GaussianModel()
    pars = gm.make_params()
    pars['center'].set(value=result.params["g0_center"].value)
    pars['amplitude'].set(value=result.params["g0_amplitude"].value)
    pars['sigma'].set(value=result.params["g0_sigma"].value)
    gm.eval(pars, x=xdata)
    ax.plot(xdata, gm.eval(pars, x=xdata), '--', color="b",label='line 1', linewidth=2)
    
    pars['center'].set(value=result.params["g1_center"].value)
    pars['amplitude'].set(value=result.params["g1_amplitude"].value)
    pars['sigma'].set(value=result.params["g1_sigma"].value) 
    gm.eval(pars, x=xdata)    
    ax.plot(xdata, gm.eval(pars, x=xdata), '--', color="r",label='line 2', linewidth=2)


def get_proj_distance( proj_pts, iq_data ):
    """
    proj_pts with shape (2,2)
    shape[0] is point idx
    shape[1] is IQ
    iq_data with shape (2,N)
    """
    p0 = proj_pts[0]
    p1 = proj_pts[1]
    ref_point = (p0+p1)/2
    shifted_iq = iq_data.transpose()-ref_point
    v_01 = p1 -p0 

    v_01_dis = np.sqrt( v_01[0]**2 +v_01[1]**2 )

    shifted_iq = shifted_iq.transpose()
    v_01 = np.array([v_01])
    projectedDistance = v_01@shifted_iq/v_01_dis
    return projectedDistance[0]

if __name__=='__main__':
    file_path = r'C:\Users\shiau\OneDrive\Desktop\multiplexRO interference'
    file_name = "IQ_Blobs_q1_2_3_4_5_all.npz"

    data = np.load(f"{file_path}\\{file_name}")# , allow_pickle=True)["arr_0"].item()
    for r in data.keys():
        
        gmm_model = train_model(data[r]*1000)
        fig = plt.figure(constrained_layout=True)
        create_img(data[r]*1000, gmm_model, fig, r)
    plt.show()