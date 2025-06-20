# %%
import numpy as np
from matplotlib import pyplot as plt
import scipy.io
from scipy import signal
from scipy.optimize import curve_fit
from ipywidgets import interact, widgets, fixed


duch_roll_data = scipy.io.loadmat("data_2024/Dutch-Roll.mat")
phugoid_data = scipy.io.loadmat("data_2024/Phugoid.mat")
roll_rub_data = scipy.io.loadmat("data_2024/Roll-Subs.mat")
spiral_data =  scipy.io.loadmat("data_2024/Spiral.mat")
SPO_data =  scipy.io.loadmat("data_2024/SPPO.mat")

duch_roll_data['Name'] = "Duch Roll"
phugoid_data['Name'] = "Phugoid"
roll_rub_data['Name'] = "Roll Subsidance"
spiral_data['Name'] = "Spiral"
SPO_data['Name'] = "SPO"

print_better = lambda stuff: print("\n".join(list(map(str, stuff))), "\n")

def secondOrd(t0,B,peak,period,zeta,t, phi=0):
    # t0,B,peak,period,zeta,t
    
    psim = np.atan(np.sqrt(1-zeta**2)/zeta)
    A = (peak - B)/(np.sqrt(1-zeta**2) * np.exp(-zeta*psim/np.sqrt(1-zeta**2)))

    omd = 2*np.pi/period
    omn = omd/np.sqrt(1-zeta**2)

    return A * np.exp(-zeta*omn*(t-t0)) * np.sin(omd*(t-t0)+phi)  +  B

# def secondOrd(t0,B,A,period,zeta,t):
#     psim = np.atan(np.sqrt(1-zeta**2)/zeta)
#     # A = (peak - B)/(np.sqrt(1-zeta**2) * np.exp(-zeta*psim/np.sqrt(1-zeta**2)))

#     omd = 2*np.pi/period
#     omn = omd/np.sqrt(1-zeta**2)

#     return A * np.exp(-zeta*omn*(t-t0)) * np.sin(omd*(t-t0))  +  B


def gauss_smooth(data, window_size=5, sigma=1):
    x = np.linspace(-1, 1, window_size)
    # print(data.shape, x.shape)
    kernel = 1/(sigma*np.sqrt(2*np.pi))*np.exp(-0.5*(x/sigma)**2)
    kernel /= np.sum(kernel)
    padded_data = np.pad(data, window_size, "edge")
    # print(np.convolve(kernel, data[:, 0],  mode="same"))
    return np.convolve(kernel, padded_data,  mode="same")[window_size:-window_size]

def plot_all(data, smooth=True, exclude=None, start=0, end=None):
    t_true = data['Time']
    if end is None : end = t_true[-1]
    try:
        start_i = np.where(start <= t_true)[0][0]
        end_i = np.where(t_true <= end)[0][-1]
    except IndexError:
        print("Invalid start or end")
        return 
    t = t_true[start_i:end_i].flatten()

    for key in data:
        if key not in  ["Time", "Name", '__header__', '__version__', '__globals__', "EAS", exclude]:
            norm_data = (data[key]-data[key][0])/np.max(np.abs(gauss_smooth((data[key]-data[key][0])[:, 0], window_size=30)))
            norm_data = norm_data[start_i:end_i].flatten()
            if smooth:
                plt.plot(t, gauss_smooth(norm_data, window_size=30), label=f"smooth {key}")
            else:
                plt.plot(t, norm_data, label=key)
    plt.title(data["Name"])
    plt.legend()
    plt.ylabel("Normalised vals")
    plt.xlabel("time (s)")
    plt.show()

def plot_one(data, key, show=True, smooth=True, start=0, end=None, norm=True):
    t_true = data['Time']
    if end is None : end = t_true[-1]
    try:
        start_i = np.where(start <= t_true)[0][0]
        end_i = np.where(t_true <= end)[0][-1]
    except IndexError:
        print("Invalid start or end")
        return 
    
    if norm:
        norm_data = (data[key]-data[key][0])/np.max(np.abs(data[key]-data[key][0]))
    else:
        norm_data = data[key][:]
        
    t = t_true[start_i:end_i].flatten()
    norm_data = norm_data[start_i:end_i].flatten()

    plt.plot(t, norm_data, label=key)
    if smooth : plt.plot(t, gauss_smooth(norm_data, window_size=30), label=f"smooth {key}")
    if show:
        plt.title(key)
        plt.ylabel("Normalised val")
        plt.xlabel("time (s)")
        plt.legend()
        plt.show()

def plot_guess(data, key, p0, phi=0, start=0, end=None, show=True):
    t_true = data['Time']
    if end is None : end = t_true[-1]
    try:
        start_i = np.where(start <= t_true)[0][0]
        end_i = np.where(t_true <= end)[0][-1]
    except IndexError:
        print("Invalid start or end")
        return 
    t = t_true[start_i:end_i].flatten()

    plt.plot(t, secondOrd(start,*p0,t, phi=phi), label="simulated")
    if show:
        plt.legend()
        plt.show()

def freq_analysis(data, key, start=0, end=None, window_size=30, plot = True, show=True, debug=False, x_range=None, ignore_first=0):
    t_true = data['Time']
    if end is None : end = t_true[-1]
    try:
        start_i = np.where(start <= t_true)[0][0]
        end_i = np.where(t_true <= end)[0][-1]
    except IndexError:
        print("Invalid start or end")
        return 
    t = t_true[start_i:end_i].flatten()
    data_n = data[key][start_i:end_i].flatten()

    T= np.mean(np.diff(t, axis=0))
    fft_domain = np.fft.fft(gauss_smooth(data_n, window_size=window_size))
    freq = np.fft.fftfreq(len(fft_domain), d=T)
    n = ignore_first
    L = len(fft_domain)
    freq_out = np.mean(abs(freq[1+signal.find_peaks(abs(fft_domain[1:]), height=np.max(abs(fft_domain[1:])))[0]]))

    if debug:
        plt.plot(freq[1:], abs(fft_domain[1:]))
        if x_range is not None: plt.xlim(-x_range, x_range)
        if show: plt.show()

    if plot:
        plot_one(data, key, show=show, start=14, norm=False)

    return freq_out

def cost_func(sim_data, real_data):
    pass

def func(t, B, A, freq, phi, one_over_tau):
    return B + A*np.exp(-t*one_over_tau)*np.sin(2*np.pi*freq*t + phi)


# key = "EAS"
# gauss_smooth(duch_roll_data[key]/duch_roll_data[key][0])

# freq_analysis(duch_roll_data, "Rollang")
plot_all(duch_roll_data, start=14, end=22)


# %%
print("Dutch-Roll")
print_better(scipy.io.whosmat("data_2024/Dutch-Roll.mat"))

print("Phugoid")
print_better(scipy.io.whosmat("data_2024/Phugoid.mat"))

print("Roll-Subs")
print_better(scipy.io.whosmat("data_2024/Roll-Subs.mat"))

print("Spiral")
print_better(scipy.io.whosmat("data_2024/Spiral.mat"))

print("SPPO")
print_better(scipy.io.whosmat("data_2024/SPPO.mat"))

# %%
plot_all(duch_roll_data, start = 14)
plot_one(duch_roll_data, 'Rudder')

start = 14

def guess(data, key, start=0, end=None, debug=False, x_range=0.5):
    t_true = data['Time']
    if end is None : end = t_true[-1]
    try:
        start_i = np.where(start <= t_true)[0][0]
        end_i = np.where(t_true <= end)[0][-1]
    except IndexError:
        print("Invalid start or end")
        return 
    t = t_true[start_i:end_i].flatten()

    freq_out = freq_analysis(data, key, start=start,end=end, show=debug, debug=debug, plot=debug, x_range=x_range)
    if debug: print(freq_out)

    data_n = data[key][start_i:end_i].flatten()

    B = np.mean(gauss_smooth(data_n)) #/(t[-1]-14)*(t[1]-t[0])
    if debug: print(B)

    # B-=s1.3

    peaks = signal.find_peaks(gauss_smooth(data_n, window_size=30))[0]
    peaks_neg = signal.find_peaks(-gauss_smooth(data_n, window_size=30))[0]

    valid_peak = peaks
    ng_valid_peak=peaks_neg
    # for i in range(peaks.shape[0]):
    #     if start <= t[peaks_neg[i]] <= end: valid_peak.append(peaks[i])

    # ng_valid_peak = []
    # for i in range(peaks_neg.shape[0]):
    #     if start <= t[peaks_neg[i]] <= end: ng_valid_peak.append(peaks_neg[i])
    # valid_peak.sort()
    # valid_peak = np.array(valid_peak[:])
    # print(duch_roll_data["Rollang"][valid_peak])
    # plt.yscale("log")
    # print(t.flatten().shape)
    coef = np.polyfit(t[valid_peak].flatten(), np.log(data[key][valid_peak].flatten()), deg=1)
    func = np.poly1d(coef)

    p0 = [B, data[key][min(valid_peak[0], abs(ng_valid_peak[0]))], 1/freq_out, -coef[0]/(freq_out*2*np.pi)]

    if debug:
        plt.plot(t, secondOrd(start,*p0,t), label="simulated")
        plt.legend()
        plt.show()

    if debug:
        plt.plot(t[valid_peak], np.log(data[key][valid_peak]), 'x')
        plt.plot(t, func(t))
        plt.show()


    print(f"peak={coef[1]}, 1/tau = {coef[0]}")
    print(f"first peak = {data[key][valid_peak[0]]}")
    return p0

def func(t, B, A, freq, phi, one_over_tau):
    return B + A*np.exp(-t*one_over_tau)*np.sin(2*np.pi*freq*t + phi)

# B, A, freq, phi, 1/tau
# p0 = np.array([B, coef[1], freq_out, 0, -coef[0]])
# tol = 1e-3
# popt, pcov = curve_fit(func, t[t>14], duch_roll_data["Rollang"][t>14], p0, bounds=[[-np.inf, 0, freq_out-tol, -np.pi, 0], [np.inf, 100, freq_out+tol, np.pi, np.inf]])
# popt
# secondOrd(14,-0.3361946354689098-1.3,7.716,1/freq_out,0.4124759024598722//(freq_out*2*np.pi,t[t>14])




# %%
p0 = list(guess(duch_roll_data, "Rollang", start=14, debug=False))

def trial_and_error(B, peak, period, zeta):
    p= [B, peak, period, zeta]
    plot_one(duch_roll_data, "Rollang", show=False, start=14, norm=False)
    plot_guess(duch_roll_data, "Rollang", p, start=14)

p0[0]=-1.71
p0[1]=5.8
p0[3]=0.2


B = widgets.FloatSlider(value=p0[0], min=-5, max=5, step=0.01)
peak = widgets.FloatSlider(value=p0[1], min=0, max=10, step=0.1)
period = widgets.FloatSlider(value=p0[2])
zeta = widgets.FloatSlider(value=p0[3], min=0, max=1, step=0.01)
print(p0[3])

interact(trial_and_error, B=B, peak=peak, period=period, zeta=zeta)

# %%
plot_all(phugoid_data, start=20)
plot_one(phugoid_data, "Elevator")

t=phugoid_data["Time"]
freq_out = freq_analysis(phugoid_data, "Nz", start=20, debug=True, x_range=0.1)
print(freq_out)
plt.plot(t[t>20], -secondOrd(20,-0.3,0.6,1/freq_out,0.2,t[t>20]), label="simulated")

plt.legend()



freq_analysis(phugoid_data, "Ptchang", start=20)

# %%
data = phugoid_data
key = "Ptchang"
start = 20

p0 = list(guess(data, key, start=start, debug=False))

def trial_and_error(B, peak, period, zeta, phi):
    p= [B, peak, period, zeta]
    plot_one(data, key, show=False, start=start, norm=False)
    plot_guess(data, key, p,phi=phi, start=start)

p0[0]=1.43
p0[1]=19.40
# p0[2]=51.2
p0[3]=0.07
phi=-0.05

B = widgets.FloatSlider(value=p0[0], min=-5, max=5, step=0.01)
peak = widgets.FloatSlider(value=-p0[1], min=-100, max=0, step=0.1)
period = widgets.FloatSlider(value=p0[2], step=0.001)
zeta = widgets.FloatSlider(value=p0[3], min=0, max=1, step=0.01)
phi=widgets.FloatSlider(value=phi, min=-np.pi, max=np.pi, step=0.01)
print(p0[2])

interact(trial_and_error, B=B, peak=peak, period=period, zeta=zeta, phi=phi)

# %%
data = phugoid_data
key = "Nz"
start = 20

p0 = list(guess(data, key, start=start, debug=False))

def trial_and_error(B, peak, period, zeta, phi):
    p= [B, peak, period, zeta]
    plot_one(data, key, show=False, start=start, norm=False)
    plot_guess(data, key, p,phi=phi, start=start)

p0[0]=1.02
p0[1]=0.56
# # p0[2]=51.2
p0[3]=0.06

phi=1.52

B = widgets.FloatSlider(value=p0[0], min=-5, max=5, step=0.01)
peak = widgets.FloatSlider(value=p0[1], min=0, max=2, step=0.01)
period = widgets.FloatSlider(value=p0[2], step=0.001)
zeta = widgets.FloatSlider(value=p0[3], min=0, max=1, step=0.01)
phi=widgets.FloatSlider(value=phi, min=-2*np.pi, max=2*np.pi, step=0.01)
print(p0[2])

interact(trial_and_error, B=B, peak=peak, period=period, zeta=zeta, phi=phi)

# %%
plot_all(SPO_data)
plot_one(SPO_data, "Elevator")

# %%
data = SPO_data
key = "Nz"
start = 9

# p0 = list(guess(data, key, start=start, debug=True, x_range=1))

def trial_and_error(B, peak, period, zeta, phi):
    p= [B, peak, period, zeta]
    plot_one(data, key, show=False, start=start, norm=False)
    plot_guess(data, key, p,phi=phi, start=start)

p0[0]=1.04
p0[1]=1.11
p0[2]=3.73
p0[3]=0.17

phi=1.3

B = widgets.FloatSlider(value=p0[0], min=1, max=1.2, step=0.01)
peak = widgets.FloatSlider(value=p0[1], min=1, max=1.2, step=0.01)
period = widgets.FloatSlider(value=p0[2], min=3.6, max=4, step=0.001)
zeta = widgets.FloatSlider(value=p0[3], min=0, max=1, step=0.01)
phi=widgets.FloatSlider(value=phi, min=-2*np.pi, max=2*np.pi, step=0.01)
print(p0[2])

interact(trial_and_error, B=B, peak=peak, period=period, zeta=zeta, phi=phi)

# %%
data = SPO_data
key = "Alpha"
start = 9

# p0 = list(guess(data, key, start=start, debug=True, x_range=1))

def trial_and_error(B, peak, period, zeta, phi):
    p= [B, peak, period, zeta]
    plot_one(data, key, show=False, start=start, norm=False)
    plot_guess(data, key, p,phi=phi, start=start)

p0[0]=3.46
p0[1]=4.03
p0[2]=3.89
p0[3]=0.41

phi=0.89

B = widgets.FloatSlider(value=p0[0], min=3, max=4.5, step=0.01)
peak = widgets.FloatSlider(value=p0[1], min=3, max=4.5, step=0.01)
period = widgets.FloatSlider(value=p0[2], min=0, max=10, step=0.001)
zeta = widgets.FloatSlider(value=p0[3], min=0, max=1, step=0.01)
phi=widgets.FloatSlider(value=phi, min=-2*np.pi, max=2*np.pi, step=0.01)
print(p0[2])

interact(trial_and_error, B=B, peak=peak, period=period, zeta=zeta, phi=phi)

# %%


# %%
plot_all(spiral_data, start=0)
plot_one(spiral_data, "Aileron", start=0, end=30)
plot_one(spiral_data, "Rollrt", start=0, end=30)
freq_out = freq_analysis(spiral_data, "Rollrt", start=0, end=30, debug=True, x_range=0.4)

data = spiral_data
key = "Rollrt"
start = 0
end = 30

# p0 = list(guess(data, key, start=start, debug=True, x_range=1))

def trial_and_error(B, peak, period, zeta, phi):
    p= [B, peak, period, zeta]
    plot_one(data, key, show=False, start=start, norm=False)
    plot_guess(data, key, p,phi=phi, start=start)

p0[0]=1
p0[1]=1
p0[2]=1/0.3
p0[3]=0.05

phi=0

B = widgets.FloatSlider(value=p0[0], min=-1, max=2, step=0.01)
peak = widgets.FloatSlider(value=p0[1], min=-1, max=2, step=0.01)
period = widgets.FloatSlider(value=p0[2], min=0, max=10, step=0.001)
zeta = widgets.FloatSlider(value=p0[3], min=0, max=1, step=0.01)
phi=widgets.FloatSlider(value=phi, min=-2*np.pi, max=2*np.pi, step=0.01)
print(p0[2])

interact(trial_and_error, B=B, peak=peak, period=period, zeta=zeta, phi=phi)

# %%
plot_all(roll_rub_data)
plot_one(roll_rub_data, "Aileron")

# %%



