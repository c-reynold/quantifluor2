import numpy as np
import pandas as pd
from pathlib import Path
import re
import json
from ops_client.client import OpsClient
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import tkinter as tk
from tkinter import ttk, filedialog, messagebox


def query_wlo_extraction_plate(barcode):

    '''Given a barcode, queries wet lab ops for that barcode and returns a json like dict'''

    #Contact Xavier Atache for username and password for WLO query. Potentially this needs to be rewritten to allow user to enter their own username and password.
    user = ""
    password = ""
    ops_url = "https://wlo-web-biopharma-stg1.apps.dev.9qzc.p1.openshiftapps.com"
    details_url = "/rest/labware/{consumable}/content/?by=barcode&all_samples=TRUE"

    labware_barcode = barcode
    client = OpsClient(ops_url, user, password)
    data = client.rest_client.get(details_url.format(consumable=labware_barcode)).json()

    return data

def _build_sigma_y2(y: np.ndarray):
    """
    Build sigma for scipy curve_fit.

    curve_fit uses sigma as standard deviations, and minimizes:
        sum(((y - f(x)) / sigma) ** 2)

    So:
    - weighting = 1 / y^2  means sigma should be proportional to y

    For 1 / y^2 weighting, we use:
        sigma = abs(y)

    To avoid division-by-zero or tiny unstable sigma values, a floor is applied.
    """
    #We first take the absolute values so that we can find an appropriate non-0 floor.
    y_abs = np.abs(y)

    #Now we are going to look at the y values to define the lowest y-value we can have.
    #For this, we can either choose 1e-12 as a pretty small number, but if y values can be less than that, then we don't want to set every y-value that is less to the floor, because then it just flattens. But we don't want them to be too absurdly small, otherwise curve fitting will only fit the bottom point.
    #Instead, we will set the floor to be the maximum y-value*1e-12, or 1e-6, whichever is larger.
    #In other words, we want to make the floor very small compared with the biggest observed signal, but never smaller than 1e-12.
    floor = max(np.nanmax(y_abs) * 1e-12, 1e-6)


    #Here we are replacing any values lower than the floor with the floor value
    y_safe = np.maximum(y_abs, floor)

    #Now because curve_fit minimizes this function:
    #sum(((y-f(x))/sigma))**2, to get y^2 weighting, we just need sigma to be equal to y, because (y-f(x)/sigma)**2=((y-f(x))**2)/(sigma**2), so if sigma=y, then we get 1/y^2 weighting
    return y_safe

def four_param_logistic(x, bottom, top, ec50, hill_slope):
    """
    4-parameter logistic (4PL) model.
    """
    x = np.asarray(x, dtype=float)
    return bottom + (top - bottom) / (1 + (x / ec50) ** (-hill_slope))

def four_pl_weighted_curve_fit(
    x,
    y,
    model_func = four_param_logistic,
    p0=[100, 10000, 2.0, 1.2],
    bounds=(-np.inf, np.inf),
    absolute_sigma: bool = False,
    n_curve_points: int = 500
):
    """
    Fit a custom model to x/y data with optional weighting and residual plots.

    Parameters
    ----------
    model_func : callable
        Function of the form f(x, *params).
    x : array-like
        X values for fit.
    y : array-like
        Y values for fit.
    p0 : list/tuple/array, optional
        Initial parameter guesses. Should leave empty unill you have 
    weighting : str or None
        Weighting scheme. Options:
            None       -> unweighted
            "1/y^2"    -> weights proportional to 1 / y^2
            "1/y"      -> weights proportional to 1 / y
    bounds : 2-tuple
        Lower and upper parameter bounds for curve_fit.
    absolute_sigma : bool
        Passed to scipy.optimize.curve_fit.
    make_plot : bool
        Whether to generate fit and residual plots.
    n_curve_points : int
        Number of points used to draw the smooth fitted curve.
    title : str
        Plot title.

    Returns
    -------
    results : dict
        Dictionary containing:
            - popt: fitted parameter values
            - pcov: covariance matrix
            - perr: parameter standard errors
            - y_fit: fitted y values at original x points
            - residuals: y - y_fit
            - sigma: sigma array used by curve_fit
    """
  

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if x.shape != y.shape:
        raise ValueError("Number of x and y fit data coordinates don't match.")

    # Remove NaN/inf rows
    valid = np.isfinite(x) & np.isfinite(y)
    if False in valid:
        print("The following standard(s) have a non-finite value:")
        print(f"x/y pairs:{zip(x[~valid],y[~valid])}")
        raise ValueError("One or more standard values is not a finite number")


    sigma = _build_sigma_y2(y=y)

    popt, pcov = curve_fit(
        f=model_func,
        xdata=x,
        ydata=y,
        p0=p0,
        sigma=sigma,
        absolute_sigma=absolute_sigma,
        bounds=bounds,
        maxfev=10000,
    )

    y_fit = model_func(x, *popt)
    residuals = y - y_fit
    weighted_residuals = residuals / sigma

    perr = np.sqrt(np.diag(pcov))

  
    _plot_fit_and_residuals(
        model_func=model_func,
        x=x,
        y=y,
        popt=popt,
        residuals=weighted_residuals,
        sigma=sigma,
        n_curve_points=n_curve_points,
        title="Standard Curve",
    )

    return {
        "popt": popt,
        "pcov": pcov,
        "perr": perr,
        "y_fit": y_fit,
        "residuals": residuals,
        "weighted_residuals":weighted_residuals,
        "sigma": sigma,
    }

def _plot_fit_and_residuals(
    model_func,
    x,
    y,
    popt,
    residuals,
    sigma,
    n_curve_points,
    title,
):
    """
    Plot data + fit and residuals.
    """
    positive_x = x[x > 0]

    x_curve = np.geomspace(
        positive_x.min(),
        positive_x.max(),
        n_curve_points
)
    y_curve = model_func(x_curve, *popt)

    fig, axes = plt.subplots(
        2, 1, figsize=(8, 8), sharex=True,
        gridspec_kw={"height_ratios": [3, 1]}
    )

    ax_fit = axes[0]
    ax_res = axes[1]

    if sigma is not None:
        ax_fit.errorbar(
            x, y, fmt="o", capsize=3, label="Data"
        )
    else:
        ax_fit.plot(x, y, "o", label="Data")

    ax_fit.plot(x_curve, y_curve, "-", label="Fit")
    ax_fit.set_ylabel("Relative Flouresence Units")
    ax_fit.set_title(title)
    ax_fit.legend()
    ax_fit.grid(True, alpha=0.3)

    ax_res.axhline(0, linestyle="--")
    ax_res.plot(x, residuals, "o")
    ax_res.set_xlabel("Nucleic Acid Mass (ng)")
    ax_res.set_ylabel("Weighted Residual")
    ax_res.grid(True, alpha=0.3)

    
    ax_fit.set_xscale("log")
    ax_fit.set_yscale("log")

    ax_res.set_xscale("log")

    #right now, the plot is just shown. It can be saved as a PNG, but not sure exactly how to do that in the moment with the way that this function is wrapped in another function in the fitting.py script.
    plt.tight_layout()
    plt.show()

def parse_quant_file(filepath:str,matrix_barcode:str):

    #Creates empty dict to append results to
    results = {}

    #Here is our data file that we will be reading. This comes from the tecan as a strangely encoded .asc file. Use latin1 for encoding when reading the file.
    raw_data_path=Path(filepath)    

    #The tecan gives us PC1 and PC2 for both the replicates of postitive controls. This doesn't match the {sample}{replicate}_{number} format of the other samples, so we need to change the positive control to match this. For this, we use a converter
    replacement_ids = {
        "PC1": "PC1_1",
        "PC2": "PC2_1",
        "BL1": "BL1_1",
        "BL2": "BL2_1"
        }

    #Here we are loading in the dataframe but we are just loading the raw data into the dataframe. We are converting the PC1 and PC2 values using the replacement IDs above and a lambda function
    df = pd.read_csv(
    raw_data_path,
    sep="\t",
    header=None,
    names=["sample_id", "well", "raw_rfu"],
    usecols=[0, 1, 2],
    index_col=False,
    nrows=96,
    encoding="latin1",
    converters={
        "sample_id": lambda value: replacement_ids.get(value, value)
        }
    )
    
    # Confirm that the file contains the expected 96 wells
    if len(df) != 96:
        raise ValueError(f"Expected 96 rows, but found {len(df)}")


    #Now we need to add the original well columns to this dataframe. These are the well locations from the original matrix plate. These will be used later to merge sample barcodes back onto the quant data
    #create an empty dict to use for mapping 

    rows = list("ABCDEFGH")
    source_columns = range(1, 11)

    well_map = {}

    for source_matrix_column in source_columns:
        quant_plate_column = source_matrix_column  + 2
        print(f"source column:{quant_plate_column}")

        for row in rows:

            if source_matrix_column < 10:
                matrix_well = f"{row}0{source_matrix_column}"

            else:
                matrix_well = f"{row}{source_matrix_column}"

            well_map[f"{row}{quant_plate_column}"] = matrix_well

    #We are deleting H11 and H12 because that is the position of the controls and we don't want them to erroneously get assigned a barcode
    del well_map["H12"]

    #Now we will apply the mapping dict to the dataframe so we can accuratly assign tube barcodes to each sample in a later step
    df["original_matrix_well"]=df["well"].map(well_map)

     
    #Now we need to clean up the sample id column so that we can get a unique pivot key. For this, we want {sample_type}{sample_number}, and we want replicate in another column so we can have replicate information in wide format.

    #Now lets make sure that all of the original matrix wells for each sample replicate match. If not, there is a problem.



    #Last thing we need to do to this dataframe is append tube barcodes. For this, we need to login to WLO to pull the corresponding matrix rack tube IDs
    wlo_query=query_wlo_extraction_plate(matrix_barcode)


    matrix_well_list=df["original_matrix_well"].dropna().unique()

    matrix_well_map={}

    for well in matrix_well_list:
        matrix_well_map[f"{well}"]=wlo_query["result"][f"{well}"]["tube_barcode"]

    #Now lets apply the map to get tube barcodes onto the dataframe

    df["matrix_tube_barcode"]=df["original_matrix_well"].map(matrix_well_map)

    #Now we will export add the result to the return dict

    results["data"]=df

    
    #Now lets grab the other metadata information from this file
    # errors="replace" prevents odd characters such as � from crashing the read. The tecan has weird text encoding when it makes files, not sure why
    text = raw_data_path.read_text(encoding="latin1", errors="replace")

    #These are regex strings for pulling metadata information from the file above.
    patterns = {
        "quant_plate_barcode": r"^Scan Barcode:\s*(.+)$",
        "instrument_serial": r"^Instrument serial number:\s*(.+)$",
        "workspace_file": r"^(.+\.wsp)$",
        "optimal_gain": r"^Value of optimal gain:\s*(\d+)$",
        "date": r"^Date:\s*(\d{4}-\d{2}-\d{2})",
        "time": r"Time:\s*(\d{2}:\d{2}:\d{2})",
        "method": r"^(Quantifluor[^\r\n]*)$",
    }


    #HEre we are using the regex patterns above to look for metadata in the result file, then we are appending the data to the results dict so it can be accessed in the future

    for name, pattern in patterns.items():
        match = re.search(pattern, text, flags=re.MULTILINE)

        if match:
            results[name] = match.group(1)
        else:
            results[name] = None


    return results



def four_pl_interpolate_concentration(y,bottom,top,ec50,hill_slope):
    
    """
    Given a flouresence value, and fit parameters, calculate the interpolated concentration for a 4PL fit
    """
    if y is pd.NA:
        return pd.NA

    if y > top:
        return ">top"
    
    if y < bottom:
        return 0

    conc=ec50*((y-bottom)/(top-y)) ** (1/hill_slope)

    return conc 

def determine_concentration_status(row):

    interpolated_mass = row["interpolated_mass (ng)"]

    reasons = []

    if pd.isna(interpolated_mass):
        reasons.append("Interpolated Mass is missing. Check raw RFU output for non-int value, or an RFU value much larger than highest standard RFU")
    elif interpolated_mass == "<bottom":
        reasons.append("replicate 1 below fitted bottom")
    elif interpolated_mass == ">top":
        reasons.append("replicate 1 above fitted top")


    if len(reasons) > 0:
        return "; ".join(reasons)

    final_concentration = row["final_concentration (ng/ul)"]
    

    if pd.isna(final_concentration):
        return "final concentration could not be calculated. Check raw RFU output for non-int value, or an RFU value much larger than highest standard RFU"

    return "reportable"

def interpolate_raw_data(data, sample_type:str):
    
    eight_std_sample_types=["cfdna","gdna","ffpe_dna"]
    seven_std_sample_types=["ffpe_rna"]

    if sample_type not in eight_std_sample_types+seven_std_sample_types:
        raise ValueError(f"Supplied sample_type {sample_type} is not in the approved list of {eight_std_sample_types+seven_std_sample_types}")


    #For DNA, we use all 8 points to fit the data. For RNA, there is not a difference between std7 and the blank, so we drop std7. This means also that the detection limit for RNA is higher. The number of ng in each well was optimized to be the same between DNA and RNA
    #x_std_mapping tells you the number of ng of nucleic acid for each standard.


    x_std_mapping={}

    num_standard_replicates=2
    for i in range(1,num_standard_replicates+1):
        temp_map={
            f"A{i}":200,
            f"B{i}":50,
            f"C{i}":12.5,
            f"D{i}":3.125,
            f"E{i}":0.78125,
            f"F{i}":0.1953125,
            f"G{i}":0.048828125,
            f"H{i}":0
        }
        x_std_mapping.update(temp_map)

    #Here we remove the 7th standard for RNA curves
    if sample_type in seven_std_sample_types:
       del x_std_mapping["G1"]
       del x_std_mapping["G2"]

    #Here we use the index values to create a new std_mass column, so we can ensure that we have paried x and y values passed correctly to the fitting function
    data["std_mass (ng)"]=data["well"].map(x_std_mapping)

    #This filters the data dataframe for just the standards by checking to see if columns are NA. They should be NA if the index is not matched to items in x_std_mapping. It also sorts them for fun.
    standards_df=data[data["std_mass (ng)"].notna()].sort_values("std_mass (ng)", ascending=False)


    #Here we are taking the mean and stdev of each standard pair 
    grouped_standards = (
        standards_df
        .groupby("std_mass (ng)", as_index=False)
        .agg(
            mean_rfu=("raw_rfu", "mean"),
            stdev_rfu=("raw_rfu", "std")
        )
    )

    grouped_standards["percent_cv"] = (grouped_standards["stdev_rfu"] / grouped_standards["mean_rfu"])*100

    #We are adding the mean, stdev, and CV to the "data" dataframe for future reference or when tech's are having issues, but the main std fitting will continue with the grouped_standards df
    data=data.merge(grouped_standards,how="left", on="std_mass (ng)")

    #Here we actually are pulling our x and y values to pass them to the fitting function
    x=grouped_standards["std_mass (ng)"]
    y=grouped_standards["mean_rfu"]

    #Here is were we actually fit, the output of which is the fit parameters and the graph.
    fit_results=four_pl_weighted_curve_fit(x=x,y=y)


    #print results to std out. Probably not needed, but nice to have
    print("Fitted parameters:")
    for name, value, err in zip(
        ["bottom", "top", "ec50", "hill_slope"],
        fit_results["popt"],
        fit_results["perr"],
    ):
        print(f"{name:>12}: {value:.6g} ± {err:.3g}")


    #here we are pulling the fit parameters out so we aren't writing fit_results["popt"][int] all over the place. Makes it easier to read and debug and work with
    bottom_fit=fit_results["popt"][0]
    top_fit=fit_results["popt"][1]
    ec50_fit=fit_results["popt"][2]
    hill_slope_fit=fit_results["popt"][3]

    #Creates a new column with the interpolation function to add the values. We do this for both replicate 1 and replicate 2 seperately, then average and CV them later rather than fitting the average RFU. Probably no right or wrong way here, should be pretty closely mathematically identical.
    # These retain the original numeric values or strings such as
    # If the values are above or below the top_fit or bottom_fit values, it will return a text string indicating as such

    data["interpolated_mass (ng)"] = data["raw_rfu"].map( 
        lambda y: four_pl_interpolate_concentration( 
            y,
            bottom_fit, 
            top_fit, 
            ec50_fit, 
            hill_slope_fit,
            ) 
        )


    #This dictionary defines the dilution factor used for each of these processes in the dilution plate. For example, ffpe_dna has a 1:50 dilution prior to measuring (2 ul into 100 ul), RNA is 1:25, etc
    sample_dilution_correction={

        "ffpe_dna":50,
        "ffpe_rna":25,
        "cfdna":5,
        "gdna":100
        }

    #This line of code is converting the number of ng in the well to the original concentration in the sample tube. The divide by 10 comes from 10 ul of diluted sample being put into the final quant plate containing dye
    data["final_concentration (ng/ul)"] = (data["interpolated_mass (ng)"] / 10 ) * sample_dilution_correction[sample_type]


    #This entry tells you if each concentration is reportable or not, and if its not, why not. There is one main reason that the sample would not be reportable, if its raw RFU value was higher than the mean standard RFU value. This means that the sample is too high in concentration and needs to be diluted. 
    #However, the "determine_concentration_status" function also checks for any reason that the raw RFU may be NaN, this can happen if the .asc file output from the plate reader itself returns a non-int value. An example of when this can happen is if the sample is WAY WAY WAY more concetrated than the highest standard.
  
    data["final_concentration_status"] = data.apply(determine_concentration_status, axis=1)
    


    ##################################################
    #Here, we are looking at batch level QC metrics
    #This pulls the control concentration out of "data"
    control_concentration = data.loc[data["well"] == "H12", "final_concentration (ng/ul)"].iloc[0]

    #This takes the data df and pulls just the standards from it
    standards_df=data[data["std_mass (ng)"] > 0]


    
    #Here we calculate the max standard CV
    standard_max_cv=standards_df["percent_cv"].max()

    #here we calculate the residual percentages for each standard
    standards_df["standard_residual_percent"] = (
        (
            standards_df["interpolated_mass (ng)"]
            - standards_df["std_mass (ng)"]
        )
        / standards_df["std_mass (ng)"]
    ) * 100


    # Here we find the largest deviation of residual percents in the standards
    #A good QC will probably be to 
    standard_max_residual = (
        standards_df["standard_residual_percent"]
        .abs()
        .max()
    )

    ##########################################



    #We need to restructure the data dataframe so that it can be ingested by autopipeline easily. For this, we want to transpose the dataframe so that the index becomes the well, and then we want to put the current index, which are the sample type descriptions, into a column

    data=data.reset_index(names="well_number")

    #This adds the 0 before 1 for wells with single digit column numbers
    data["well"] = data["well"].str[0] + data["well"].str[1:].str.zfill(2)
    data_long=data.set_index("well")
    data_long=data_long.sort_index()

    #Now we need to fill the NaN values with "null" so that it can be ingested by autopipeline

    data_long=data_long.fillna("null")

    #creating a results dict to output results

    results={}

    #add the fit parameters to the dataframe. Might be better to tag them onto the end instead of making new columns with the values.
    results["bottom_fit"]=bottom_fit
    results["top_fit"]=top_fit
    results["ec50_fit"]=ec50_fit
    results["hill_slope_fit"]=hill_slope_fit
    results["control_concentration"]=control_concentration
    results["standard_max_cv"]=standard_max_cv
    results["standard_max_residual"]=standard_max_residual
    
    
    #Adding results dataframe to function output
    results["quant_data"]=data_long.transpose().to_dict()

    final_result={}
    final_result["dict"]=results
    final_result["dataframe"]=data_long

    return final_result






def run_quant_program(
    input_file: str,
    output_file: str,
    qnt_plate_barcode: str,
    matrix_plate_barcode:str,
    sample_type:str
) -> None:
    """
    """

    # Example:
    # raw_data = pd.read_csv(input_file)
    # results = analyze_quant_plate(raw_data, plate_barcode)
    

    print(f"Input: {input_file}")
    print(f"Output: {output_file}")
    print(f"Quant Plate Barcode: {qnt_plate_barcode}")
    print(f"Matrix Plate Barcode: {matrix_plate_barcode}")
    print(f"Sample type: {sample_type}")

    parse_quant_output=parse_quant_file(input_file,matrix_barcode=matrix_plate_barcode)

    print("Sample parsing successful. Checking parsed quantities...")

    parsed_qnt_plate_barcode=parse_quant_output["quant_plate_barcode"]
    parsed_method=parse_quant_output["method"]

    print(f"Parsed quant plate barcode: {parsed_qnt_plate_barcode}")
    print(f"Parsed tecan method: {parsed_method}")

    parsed_qnt_plate_barcode = str(
        parse_quant_output["quant_plate_barcode"]
    ).strip()

    parsed_method = str(
        parse_quant_output["method"]
    ).strip()

    print(f"Parsed quant plate barcode: {parsed_qnt_plate_barcode!r}")
    print(f"Parsed Tecan method: {parsed_method!r}")

    method_to_sample_type = {
        "Quantifluor 2 DNA.mth": "ffpe_dna",
        "Quantifluor 2 RNA.mth": "ffpe_rna",
        "Quantifluor 2 cfDNA.mth": "cfdna",
        "Quantifluor 2 gDNA.mth": "gdna",
    }

    try:
        parsed_sample_type = method_to_sample_type[parsed_method]
    except KeyError:
        raise ValueError(
            f"Unrecognized Tecan method: {parsed_method!r}"
        )

    print(f"Sample type indicated by method: {parsed_sample_type}")

    if qnt_plate_barcode.upper().strip() != parsed_qnt_plate_barcode.upper():
        raise ValueError(
            f"Scanned quant plate barcode ({qnt_plate_barcode}) does not match "
            f"the barcode present in the uploaded file "
            f"({parsed_qnt_plate_barcode})."
        )

    if sample_type.strip().lower() != parsed_sample_type:
        raise ValueError(
            f"{parsed_method} was used to generate the quant file, which "
            f"corresponds to {parsed_sample_type}, but the user selected "
            f"{sample_type}, like an idiot."
        )

    final_result=interpolate_raw_data(data=parse_quant_output["data"],sample_type=sample_type)

    with open(output_file, "w") as json_file:
        json.dump(final_result["dict"], json_file, indent=4)

    final_result["dataframe"].to_csv(f"{output_file}.csv", index=False)


# #DON"T FORGET TO DELETE THIS BLOCK!!!
# run_quant_program(
#     input_file="development_test_inputs/RAS-LP-27Aug26-cfDNA3RNA3.asc",
#     output_file="test_output.json",
#     qnt_plate_barcode="RAS-LP-27Aug26-cfDNA3RNA3",
#     matrix_plate_barcode="QRS_DNA_000000",
#     sample_type="ffpe_dna"

# )



def select_input_file() -> None:
    selected_file = filedialog.askopenfilename(
        title="Select Quantification File",
        filetypes=[
            ("ASC Files", "*.asc")
        ],
    )

    if selected_file:
        input_path_var.set(selected_file)


def select_output_file() -> None:
    selected_file = filedialog.asksaveasfilename(
        title="Save Quantification Results",
        defaultextension=".json",
        filetypes=[
            ("JSON files", "*.json"),
            ("All files", "*.*"),
        ],
    )

    if selected_file:
        output_path_var.set(selected_file)


def submit() -> None:
    input_file = input_path_var.get().strip()
    output_file = output_path_var.get().strip()
    matrix_plate_barcode = matrix_barcode_var.get().strip()
    qnt_plate_barcode = qnt_barcode_var.get().strip()
    sample_type = sample_type_var.get().strip()

    if not input_file:
        messagebox.showerror(
            "Missing Input",
            "Select an input file."
        )
        return

    if not Path(input_file).is_file():
        messagebox.showerror(
            "Invalid Input",
            "The selected input file does not exist."
        )
        return

    if not output_file:
        messagebox.showerror(
            "Missing Output",
            "Select an output file."
        )
        return

    if not matrix_plate_barcode:
        messagebox.showerror(
            "Missing the Matrix Plate Barcode",
            "Enter the matrix plate barcode."
        )
        return
    
    if not qnt_plate_barcode:
        messagebox.showerror(
            "Missing Quant Plate Barcode",
            "Enter the quant plate barcode."
        )
        return
    
    if not sample_type:
        messagebox.showerror(
            "Missing Sample Type",
            "Please select sample type."
        )
        return

    try:
        status_var.set("Running...")
        root.update_idletasks()

        run_quant_program(
            input_file=input_file,
            output_file=output_file,
            qnt_plate_barcode=qnt_plate_barcode,
            matrix_plate_barcode=matrix_plate_barcode,
            sample_type=sample_type
        )

        status_var.set("Complete")

        messagebox.showinfo(
            "Quantification Complete",
            f"Results were saved to:\n{output_file}"
        )

    except Exception as error:
        status_var.set("Failed")

        messagebox.showerror(
            "Quantification Failed",
            str(error)
        )


root = tk.Tk()
root.title("Chris's Quantifluor 2: The Electric Boogaloo")
root.geometry("750x260")

main_frame = ttk.Frame(root, padding=20)
main_frame.pack(fill="both", expand=True)

input_path_var = tk.StringVar()
output_path_var = tk.StringVar()
matrix_barcode_var = tk.StringVar()
qnt_barcode_var = tk.StringVar()
sample_type_var = tk.StringVar()
status_var = tk.StringVar(value="Ready")


# Input file
ttk.Label(
    main_frame,
    text="Input file:"
).grid(row=0, column=0, sticky="w", pady=5)

ttk.Entry(
    main_frame,
    textvariable=input_path_var,
    width=70
).grid(row=0, column=1, sticky="ew", padx=10)

ttk.Button(
    main_frame,
    text="Browse",
    command=select_input_file
).grid(row=0, column=2)


# Output file
ttk.Label(
    main_frame,
    text="Output file:"
).grid(row=1, column=0, sticky="w", pady=5)

ttk.Entry(
    main_frame,
    textvariable=output_path_var,
    width=70
).grid(row=1, column=1, sticky="ew", padx=10)

ttk.Button(
    main_frame,
    text="Browse",
    command=select_output_file
).grid(row=1, column=2)


# Matrix Plate Barcode
ttk.Label(
    main_frame,
    text="Matrix Plate barcode:"
).grid(row=2, column=0, sticky="w", pady=5)

ttk.Entry(
    main_frame,
    textvariable=matrix_barcode_var
).grid(row=2, column=1, sticky="ew", padx=10)



#Quant Plate Barcode
ttk.Label(
    main_frame,
    text="Quant Plate Barcode:"
).grid(row=3, column=0, sticky="w", pady=5)

ttk.Entry(
    main_frame,
    textvariable=qnt_barcode_var
).grid(row=3, column=1, sticky="ew", padx=10)


#Sample dropdown
sample_type_options = [
    "ffpe_dna",
    "ffpe_rna",
    "cfdna",
    "gdna",
]

ttk.Label(
    main_frame,
    text="Sample Type:"
).grid(row=4, column=0, sticky="w", pady=5)

sample_type_combo = ttk.Combobox(
    main_frame,
    textvariable=sample_type_var,
    values=sample_type_options,
    state="readonly"
)

sample_type_combo.grid(
    row=4,
    column=1,
    sticky="ew",
    padx=10
)



# Run button
ttk.Button(
    main_frame,
    text="Run Quantification",
    command=submit
).grid(row=5, column=1, pady=20)


# Status
ttk.Label(
    main_frame,
    textvariable=status_var
).grid(row=6, column=1)


main_frame.columnconfigure(1, weight=1)

root.mainloop()

