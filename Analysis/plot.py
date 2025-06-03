import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import fiber_set as fs
import matplotlib.patches as mpatches
from scipy.interpolate import interp1d

font = {'size'   : 24}

matplotlib.rc('font', **font)

hbar_evpj = 6.626e-34 * 6.242e18  # Convert Hz to eV
hc = 1239.84193 # Convert eV to nm

def frame_spectra(fiber_set, plist=[], mlist=[]):
    wls = []
    Is = []
    Ps = []
    params = {param : [] for param in plist}
    metrics = {metric : [] for metric in mlist}
    for run in fiber_set.runs:
        wl = run.make_wavelength_scale()
        wls.append(wl)
        #NORMILAZATION IS TAKEN CARE HERE
        I = run.apply_jacob(normed = True)[-1]        
        Is.append(I)
        for param in plist:
            params[param] = np.concatenate([np.repeat(run.params[param], len(wl)), params[param]])
        for metric in mlist:
            metrics[metric] = np.concatenate([np.repeat(run.metrics[metric], len(wl)), metrics[metric]])
        #Added Phase functionality
        P = run.phases()[-1]
        Ps.append(P)
    dic = {"wl" : np.concatenate(wls),
            "I" : np.concatenate(Is),
            "P" : np.concatenate(Ps),
          }
    if plist:
        dic.update(params)
    if mlist:
        dic.update(metrics)
    return pd.DataFrame(dic)

def frame_full_spectra(fiber_run):
    wl = fiber_run.make_wavelength_scale()
    Is = fiber_run.apply_jacob(normed = True)
    dic = {"wl" : np.tile(wl, len(fiber_run)), "I" : np.concatenate(Is)}
    return pd.DataFrame(dic)


def frame_pulse(fiber_set, plist, mlist):
    t = np.tile(fiber_set.runs[0].make_time_scale(), len(fiber_set.runs))
    I = np.concatenate([np.power(np.abs(run.fields[-1]),2) for run in fiber_set.runs])
    dic = {"t": t, "I" : I}

    #Make columns given by plist and mlist.
    dic.update({param : [run.params[param] for run in fiber_set.runs]
                for param in plist})
    dic.update({metric : [run.metrics[metric] for run in fiber_set.runs]
                for metric in mlist})

    return pd.DataFrame(dic)

def frame(fiber_set, plist, mlist):
    dic = {}
    #Make columns given by plist and mlist.
    dic.update({param : [run.params[param] for run in fiber_set.runs]
                for param in plist})
    dic.update({metric : [run.metrics[metric] for run in fiber_set.runs]
                for metric in mlist})
    return pd.DataFrame(dic)


def plot_spect_grid(fiber_set, p1, p2=None, excel_data_path=None, phase=False, save_txt_path=None, data_on_top=False):
    """
        plot_spect_grid(fiber_set, p1, p2=None, excel_data_path=None, phase=False, save_txt_path=None, data_on_top=False)

    p1 (str): Each unique value of p1 becomes a separate subplot.
    p2 (str, optional): If provided, curves for each p2 value are plotted in different colors.
    excel_data_path (str, optional): Path to an Excel file containing experimental data to overlay on top of the simulation.
    phase (bool, default=False): Plots the spectral phase using a secondary y-axis.
    save_txt_path (str, optional): Saves wavelength and intensity data (between 550 and 1100 nm) to a text file.
    data_on_top (bool, default=False): Stacks all spectra on top of each other with vertical offsets.

    """
        # importing lab data and normalizing
    if excel_data_path:
        excel_data = pd.read_excel(excel_data_path, usecols=[0, 1], skiprows=6, names=["Wavelength", "Intensity"])
        excel_data["Intensity"] /= excel_data["Intensity"].max()

    custom_palette = ["#FF0000",  # Red
                      "#FFA500",  # Orange
                      "#90EE90",  # Light Green
                      "#008000",  # Green
                      "#40E0D0",  # Turquoise
                      "#0000FF",  # Blue
                      "#00008B"]  # Deep Blue
    
    df = frame_spectra(fiber_set, plist=[p1] + ([p2] if p2 else []), mlist=["l_edge", "r_edge"])
    unique_p1 = sorted(df[p1].unique()) #for legend
    if p2:
        g = sns.FacetGrid(df, row=p1, aspect=1.5, height=7, hue=p2, palette=custom_palette, row_order=unique_p1,sharex=False)
        unique_p2 = df[p2].unique()  # Get unique values of p2 for legend
        # Somtimes you have to use reversed(unique_p2) or unique_p2
        # Idk why the assinging of p2 and the colors is sometimes wrong
        color_map = {val: color for val, color in zip(reversed(unique_p2), custom_palette)}
    else:
        g = sns.FacetGrid(df, row=p1, aspect=2, height=10, hue=p1, palette=custom_palette)


    # Making my own Titles as the default flips the order
    for ax in g.axes.flat:
        ax.set_title("")
    # Manually add titles with correct Pin values in ascending order
    for i, pin_value in enumerate(reversed(unique_p1)):
        g.axes[i, 0].text(0.5, 1.05, f"{p1} = {pin_value} uJ", ha='center', va='center', transform=g.axes[i, 0].transAxes, fontsize=32)

    xmin = 0.95 * df["l_edge"].min(axis=0)
    xmax = 1.05 * df["r_edge"].max(axis=0)
    #Rounding xmin down, xmax up to the nearest multiple of 50 to make it look nicer
    xmin = 50 * (xmin // 50)
    xmax = 50 * ((xmax + 49) // 50)
    ymin = 0
    ymax = 0
    if phase:
        global_phase_min = df["P"].min()
        global_phase_max = df["P"].max()
    # Define a plotting function with an offset for each row
    def plot_with_offset(data, color, **kwargs):
        nonlocal ymax
        offset = ymax#use the current ymax as the offset
        if data_on_top==True:
            plt.plot(data["wl"], data["I"]+offset, linewidth=3, color=color)
        else:
            plt.plot(data["wl"], data["I"], linewidth=3, color=color)
        # Plotting the phase
        if phase:
            ax = plt.gca()  # Get current axis
            ax2 = ax.twinx()  # Create secondary y-axis
            ax2.plot(data["wl"], data["P"], linestyle="dashed", color="black", alpha=0.7, linewidth=2, label="Phase")  
            ax2.set_ylabel("Phase (rad)", color="black")  
            ax2.tick_params(axis="y", labelcolor="black")
            ax2.set_ylim(global_phase_min, global_phase_max)
        # Update ymax for the next row
        ymax += data["I"].max()
        #for plotting lab data ontop of simdata
        if excel_data_path:
            plt.plot(excel_data["Wavelength"], excel_data["Intensity"], linewidth=2, color="black", linestyle="--")
        plt.grid(True)
        # if you want to save x and y data as a txt file
        if save_txt_path:
            # Filter to include only data within xmin and xmax
            filtered = data[(data["wl"] >= 550) & (data["wl"] <= 1100)]
            with open(save_txt_path, "a") as f:
                current_row = data.iloc[0]
                f.write(f"# {p2} = {current_row[p2]}\n")
                np.savetxt(
                    f,
                    np.column_stack((filtered["wl"], filtered["I"])),
                    fmt="%.6f",
                    delimiter="\t",
                    header="Wavelength (nm)\tIntensity",
                    comments=""
                )
                print("Test")
    # Map the modified plotting function onto the FacetGrid
    g.map_dataframe(plot_with_offset)
    for ax in g.axes.flat:
        ax.spines['top'].set_visible(True)
        ax.spines['right'].set_visible(True)
        ax.spines['left'].set_visible(True)
        ax.spines['bottom'].set_visible(True)
    g.set(xlim=(xmin, xmax))
    g.set(xticks=range(int(xmin), int(xmax) + 1, 50))
    #g.set(ylim=(0, ymax+.1))
    g.set_xlabels("Wavelength (nm)")
    g.set_ylabels("Intensity (arb.)")
    for ax in g.axes.flat:
        ax.tick_params(axis='x', labelsize=14)
    # Create custom legend with colors corresponding to each unique p2
    if p2:
        legend_handles = [mpatches.Patch(color=color_map[val], label=str(val)) for val in unique_p2]
        legend_handles.reverse()#so that the order of the label and plot matches
        g.add_legend(handles=legend_handles,title=p2,bbox_to_anchor=(1.1,0.5))
    #g.savefig('Plot.png', dpi=300)

    if excel_data_path:
    # Add legend for Sim Data and Lab Data
        handles = [
            #plt.Line2D([0], [0], color='#FF0000', linewidth=2, label='Sim Data'),
            plt.Line2D([0], [0], color='black', linestyle='--', linewidth=2, label='Lab Data')
        ]
        g.add_legend(handles=handles, title="", loc='lower right')
    plt.show()
    return g

def plot_full_spectrum(fiber_run, save_txt_path=None):
    wl = fiber_run.make_wavelength_scale()
    Is = np.array(fiber_run.apply_jacob(normed = True))
    z = np.linspace(0, fiber_run.params["zmax"], len(fiber_run.fields))
    X, Y = np.meshgrid(z, wl)
    fig, ax = plt.subplots(1)
    ax.pcolormesh(X, Y, Is.T, figure=fig, cmap="viridis")#viridis, inferno, gist_heat,nipy_spectral
    l_edge = fiber_run.metrics["l_edge"] * 0.95
    r_edge = fiber_run.metrics["r_edge"] * 1.05
    ax.set_ylim(l_edge, r_edge)
    ax.set_xlabel("Propagation Length (m)")
    ax.set_ylabel("Wavelength (nm)")

    # Optional: Save filtered data based on the y-limits (Wavelength range)
    if save_txt_path:
        # Find the indices of wavelengths that are within the specified range
        valid_indices = (wl >= l_edge) & (wl <= r_edge)

        # Filter the wavelengths and intensities to save
        filtered_wl = wl[valid_indices]
        filtered_Is = Is[:, valid_indices]  # Keep the corresponding intensities for valid wavelengths

        # Save the filtered data
        with open(save_txt_path, "w") as f:
            f.write("# Wavelengths (rows) vs Propagation Distance (columns)\n")
            f.write("# Each row corresponds to one wavelength\n")
            f.write("# First column: Wavelength (nm), followed by intensity values at each z\n")
            
            # Save only the data within the specified range
            for i in range(len(filtered_wl)):
                row = [f"{filtered_wl[i]:.6f}"] + [f"{val:.6e}" for val in filtered_Is.T[i]]
                f.write("\t".join(row) + "\n")
    plt.show()
    return ax

def metric_heatmap(fiber_set, p1, p2, metric):
    df = frame(fiber_set, plist=[p1,p2], mlist=[metric]).pivot(index=p1, columns=p2, values=metric)
    return sns.heatmap(df)

def plot_before_final(fiber_run):
    t = fiber_run.make_time_scale()
    wl = fiber_run.make_wavelength_scale()

    fig, (ax1,ax2) = plt.subplots(1, 2, sharey=True)

    # Time Plot
    Ii = np.power(np.abs(fiber_run.fields[0]),2)
    Ii = Ii/max(Ii)
    If = np.power(np.abs(fiber_run.fields[-1]),2)
    If = If/max(If)
    ax1.plot(t, Ii)
    ax1.plot(t, If)
    tmin = -fiber_run.metrics["l_edge"] / 2 * 1.05
    tmax = fiber_run.metrics["l_edge"] / 2 * 1.05
    ax1.set_xlim(tmin,tmax)
    ax1.set_xlabel("Time (fs)")
    ax1.set_ylabel("Intensity (A.U.)")
    
    # Spectrum Plot
    spectra = fiber_run.apply_jacob(True)
    Fi = spectra[0]
    Ff = spectra[-1]
    ax2.plot(wl, Fi, label="Initial")
    ax2.plot(wl, Ff, label="Final")
    wlmin = 0.95 * fiber_run.metrics["l_edge"]
    wlmax = 1.05 * fiber_run.metrics["r_edge"]
    ax2.set_xlim(wlmin, wlmax)
    ax2.set_ylim(0,1.05)
    ax2.set_xlabel("Wavelength (nm)")
    plt.legend(bbox_to_anchor=(0.62, 0.95), loc=2, borderaxespad=0.)
    return fig

def plot_field_grid(fiber_set, p1, p2):
    for run in fiber_set.runs:
        param_value1 = run.params.get(p1, None)
        param_value2 = run.params.get(p2, None)
        run.plot_field(normed=True,p1=param_value1, p2=param_value2)
