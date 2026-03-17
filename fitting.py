import numpy as np
import pandas as pd
from curve_fit_utils import weighted_curve_fit
from curve_fit_utils import four_pl_interpolate_concentration
from pathlib import Path


def four_param_logistic(x, bottom, top, ec50, hill_slope):
    """
    4-parameter logistic (4PL) model.
    """
    x = np.asarray(x, dtype=float)
    return bottom + (top - bottom) / (1 + (x / ec50) ** (-hill_slope))



def interpolate_concentrations(file):

    #read the original CSV file and load into a dataframe
    df=pd.read_csv(file, encoding="latin1", sep="\t")

    #Removes data that is not needed for cleaner experience.
    #Uses "Blank reduction" as the tecan already does the blank analysis, it might be better long term to have blank subraction occur in the script so we have control over how it happens.
    df_small=df.loc[0:95,["Layout","Original Concentrations","Blank reduction"]]

    #creates another datafram wihtout any of the samples or other junk in it, just the standards
    standards=df_small.iloc[0:16,:]

    #takes the average of the two standards and makes a datafram from them
    averaged_standards=standards.groupby(["Layout","Original Concentrations"])["Blank reduction"].mean().reset_index()


    #real example data
    #Right now, I don't think that the 0,0 point, i.e. the blank value/0 standard value is being added to the analysis.
    #adding it or not isn't a big difference, but adding it might make the curve fit look just a titch better.
    #It is not included now because there is no "original concentration" value for the balnks, so the groupby function removes it.
    x = np.array(averaged_standards["Original Concentrations"], dtype=float)

    y = np.array(averaged_standards["Blank reduction"], dtype=float)

    # # Example fake data
    #You can use these data to test if a file is not available
    # x = np.array([0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50], dtype=float)
    # y = np.array([120, 160, 310, 520, 780, 980, 1080, 1110, 1120], dtype=float)

    # Initial guesses
    p0 = [100, 10000, 2.0, 1.2]

    #do the curve fitting
    results = weighted_curve_fit(
        model_func=four_param_logistic,
        x=x,
        y=y,
        p0=p0,
        weighting="1/y^2",
        make_plot=True,
        bounds=([0, 0, 1e-9, 0.01], [np.inf, np.inf, np.inf, 10]),
        title="4PL Fit with 1/y^2 Weighting",
    )

    #print results to std out
    print("Fitted parameters:")
    for name, value, err in zip(
        ["bottom", "top", "ec50", "hill_slope"],
        results["popt"],
        results["perr"],
    ):
        print(f"{name:>12}: {value:.6g} ± {err:.3g}")


    #define the fit parameters
    bottom_fit=results["popt"][0]
    top_fit=results["popt"][1]
    ec50_fit=results["popt"][2]
    hill_slope_fit=results["popt"][3]

    #Creates a new column with the interpolation function to add the values.
    df["final_result"] = four_pl_interpolate_concentration(df["Blank reduction"],bottom_fit,top_fit,ec50_fit,hill_slope_fit)

    #add the fit parameters to the dataframe. Might be better to tag them onto the end instead of making new columns with the values.
    df["bottom_fit"]=bottom_fit
    df["top_fit"]=top_fit
    df["ec50_fit"]=ec50_fit
    df["hill_slope_fit"]=hill_slope_fit


    return df
    
#define the input and output directories
input_directory=Path("raw_data")
output_directory=Path("analyzed_data")

for file in input_directory.iterdir():

    #get the filename
    filename=file.stem

    interpolate_concentrations(file).to_csv(output_directory / f"{filename}.analyzed.csv")