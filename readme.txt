This document goes over usage of the quantifluor2_fitting scripts:
    quantifluor2_fitting_duplicate.py
    quantifluor2_fitting_singlet.py

Packages used during development/testing:
    counsyl_autoapp_client 1.0.3
    counsyl_xmlrpc         1.2.3
    cycler                 0.12.1
    fonttools              4.62.1
    idna                   3.18
    kiwisolver             1.5.0
    matplotlib             3.10.8
    numpy                  2.4.3
    ops_client             2.7.3
    packaging              26.0
    pandas                 3.0.1
    pillow                 12.1.1
    pip                    24.0
    pyparsing              3.3.2
    python-dateutil        2.9.0.post0
    pytz                   2023.3.post1
    requests               2.27.0
    scipy                  1.17.1
    six                    1.17.0
    urllib3                1.26.20

Purpose:

The purpose of these scripts are to turn raw DNA quantificaiton data (RFU) from the Tecan Infinite into nucleic acid quants (ng/ul)
To do this, plates are prepared in the laboratory according to quantifluor instructions:
    https://mygn.atlassian.net/wiki/spaces/RBD/pages/edit-v2/3201433993

Then, this program is used to analyze the .asc output.

Currently, there are 3 important pieces of informaiton from the .asc file that are used in these scripts.

    1. Sample ID, Well, RFU:
        These are the first 96 rows and contain the RFU value for each well along with a sample id.
        In the singlet script, the sample ID doesn't matter, the wells are used instead of the sample ID. The update to the duplicate script for this has not been made yet.

    2. Line 98 (Quantifluor 2 DNA.mth) 
        shows the method that was used on the tecan to actually run the quant.
        The program checks to make sure that the correct program was used because RNA and DNA have different settings that need to be used (excitiation and emission wavelength)
    3. Barcode (Scan Barcode: singlet_test)
        This line shows the quant plate barcode that was used. The magellen method requires a usre to scan in the quant plate barcode.
        This barcode is then attached to the end of the output file to ensure that users are not analyzing the incorrect quant plate with this script.
    
    There are currently checks in this script to make sure that the user is selecting not only the correct method when analyzing, but also selects the correct quant plate so that sample quants can't get swapped.

########QUICKSTART#################

Use:

    ***************************IMPORTANT**********************************
    The script will not currently run as is. The query_wlo_extraction_plate function requires a username and password input for AP. For this password, please reach out to Xavier Atache.
    There is probably a way to have users that have access to WLO enter their own username and password, but I have not been able to get user accounts to work because it seems that they lack permissions.
    So instead, I was using a username and password provided by Xavier.
    ***********************************************************************

    To use the script, just run the script you want (singlet or duplicate) and a UI should pop up.

    The UI asks for several things:
        Input File: This should be the path to the .asc file (output from tecan) that the user wants to analyze
        Output File: THis is the path of the output file.

            *Note- known bug: If the user just types the path into the box without the file suffix, the output file will not have a suffix. (example: .json)
            If the user clicks the "browse" button, the .json suffix is automatically added.
            Additionally, I have it exported currently as a .json for autopipeline and as a .csv so that the users can look at values
                The .csv suffix for this is automatically added, and it makes weird .json.csv paths. Its just a bug I didn't have time to fix. Should be trivial to strip suffixes from filename variables and then tag on a .csv and .json to the respective exports.
        Matrix Plate Barcode: THis is the barcode of the matrix plate that samples are being taken out of to quant. This barcode allows the WLO query to append tube barcodes to the samples so that the laboratory folk can easily find/veryify tubes/problems
        Quant Plate Barcode: THis is the barcode of hte quant plate itself, with the quantifluor solution and diluted samples in it. The entry here should match the barcode on the .asc file. If not, the program will flag and error and will not analyze.
        Sample type: The sample type that is being run. The difference in the sample type is twofold:
            For DNA sample types:
                The different DNA sample types have different expected dilution ratios:
                    cfDNA 1:10
                    FFPE DNA: 1:50
                    gDNA (buffy): 1:100
            For RNA samples:
                The botom standard curve point is dropped (only 7 standards used for fitting): See confluence page above.
                The RNA dilution factor is 1:25
        
    After the settings are correctly inputed, the program will fit the standards and show a standard curve for the user. If there are visual problems with the standards, the user can use this graph to try and find out what happened.
    AFter this graph is closed, the output .json and .csv files will be output to the path the user indicated.

# Quantifluor 2 Quantification Program

## Overview

This program automates the processing of Quantifluor nucleic-acid quantification runs generated on a Tecan plate reader.

The goal of the program is to take the raw fluorescence output from the instrument and convert it into a structured, traceable set of nucleic-acid concentration results that can be reviewed by a user or consumed by downstream laboratory software.

The program performs several operations that would otherwise need to be done manually:

1. Reads the raw Tecan `.asc` output file.
2. Extracts plate and instrument metadata from the file.
3. Associates quantification-plate wells with the corresponding wells on the source matrix plate.
4. Queries Wet Lab Ops (WLO) to retrieve the tube barcode associated with each matrix-plate well.
5. Identifies the Quantifluor standards and their known nucleic-acid masses.
6. Calculates replicate statistics for the standards.
7. Fits a weighted four-parameter logistic (4PL) standard curve.
8. Uses the fitted curve to interpolate the nucleic-acid mass of unknown samples.
9. Corrects the interpolated values for assay volume and process-specific dilution.
10. Calculates standard-curve and control QC metrics.
11. Performs basic process checks to prevent mismatched plate barcodes or assay methods.
12. Exports the results as JSON and CSV files.

The application also contains a Tkinter graphical interface so that the analysis can be run without directly interacting with the Python functions.

---

# Why This Program Exists

Fluorescence-based nucleic-acid quantification requires more than simply reading the RFU value produced by the plate reader.

The fluorescence signal must be interpreted relative to a set of standards containing known amounts of nucleic acid. The relationship between nucleic-acid mass and fluorescence is not assumed to be perfectly linear across the full measurement range, so this program models the standard curve using a four-parameter logistic function.

In addition, laboratory quantification results need to retain sample identity and process context. A fluorescence value is not particularly useful by itself unless it can be associated with the correct sample tube, plate, assay method, and dilution scheme.

This program therefore combines three types of information:

* **Instrument data** from the Tecan `.asc` file
* **Sample identity information** from Wet Lab Ops
* **Quantitative analysis** from the fitted standard curve

The result is a single analysis workflow that converts raw instrument output into sample-level concentration results with associated QC and metadata.

---

# General Workflow

```text
User selects Tecan .asc file
        |
        v
Parse raw fluorescence data
        |
        v
Extract plate and instrument metadata
        |
        v
Validate plate barcode and Tecan method
        |
        v
Query Wet Lab Ops using matrix plate barcode
        |
        v
Map quant wells -> matrix wells -> tube barcodes
        |
        v
Identify Quantifluor standards
        |
        v
Calculate standard replicate mean, SD, and CV
        |
        v
Fit weighted four-parameter logistic curve
        |
        v
Interpolate sample RFU values
        |
        v
Correct for dilution and assay volume
        |
        v
Calculate QC metrics
        |
        v
Export JSON and CSV results
```

---

# Input Data

## Tecan Quantification File

The primary input is a `.asc` file produced by the Tecan plate reader.

The first 96 rows of the file are interpreted as the 96 wells of the quantification plate and are read into a Pandas DataFrame containing:

* `sample_id`
* `well`
* `raw_rfu`

The program verifies that exactly 96 rows were recovered from the file.

Some Tecan identifiers for controls do not follow the same naming convention as the other samples. These values are standardized during parsing. For example:

```text
PC1 -> PC1_1
PC2 -> PC2_1
BL1 -> BL1_1
BL2 -> BL2_1
```

---

# Metadata Extraction

In addition to the well-level fluorescence data, the program reads the full `.asc` file as text and extracts run metadata using regular expressions.

The extracted metadata includes:

* Quantification plate barcode
* Instrument serial number
* Tecan workspace file
* Optimal gain
* Run date
* Run time
* Tecan Quantifluor method

This information is useful both for traceability and for verifying that the uploaded file corresponds to the assay the user intends to analyze.

---

# Sample Identity and Wet Lab Ops Integration

The Tecan output does not provide all of the sample identity information needed for downstream analysis.

The program therefore queries Wet Lab Ops using the source matrix plate barcode.

```python
query_wlo_extraction_plate(barcode)
```

The WLO query returns the contents of the matrix plate, including the tube barcode occupying each well.

The program constructs a mapping between:

```text
Quantification plate well
        ->
Original matrix plate well
        ->
Matrix tube barcode
```

The resulting tube barcode is added to the quantification DataFrame as:

```text
matrix_tube_barcode
```

This allows a measured RFU and calculated concentration to remain associated with the physical sample from which the measurement originated.

---

# Standard Curve

## Standard Masses

The Quantifluor plate contains standards with known quantities of nucleic acid.

The DNA standard masses used by the program are:

| Standard |           Mass |
| -------- | -------------: |
| A        |         200 ng |
| B        |          50 ng |
| C        |        12.5 ng |
| D        |       3.125 ng |
| E        |     0.78125 ng |
| F        |   0.1953125 ng |
| G        | 0.048828125 ng |
| H        |           0 ng |

Two replicates are measured for each standard.

For FFPE RNA, the G standard is excluded from the curve because the low-concentration standard is not sufficiently distinguishable from the blank in the current assay configuration.

---

# Standard Replicate Statistics

Before fitting the calibration curve, replicate measurements for each standard concentration are grouped together.

For each standard, the program calculates:

* Mean RFU
* Standard deviation of RFU
* Percent coefficient of variation (%CV)

The mean RFU of the replicate pair is used as the fluorescence value for curve fitting.

The replicate statistics are also retained in the output to support troubleshooting and QC review.

---

# Four-Parameter Logistic Model

The standard curve is modeled using a four-parameter logistic, or 4PL, function:

```text
                    top - bottom
y = bottom + ------------------------------
               1 + (x / EC50)^(-Hill)
```

where:

* `x` = nucleic-acid mass
* `y` = predicted fluorescence
* `bottom` = lower asymptote
* `top` = upper asymptote
* `EC50` = mass corresponding approximately to the midpoint of the response range
* `hill_slope` = steepness of the response curve

The implementation is:

```python
bottom + (top - bottom) / (1 + (x / ec50) ** (-hill_slope))
```

The fitted parameters are retained as part of the run output.

---

# Weighted Curve Fitting

The standard curve is fit using `scipy.optimize.curve_fit`.

A `1/y^2` weighting scheme is used.

`curve_fit` minimizes:

```text
sum(((observed - predicted) / sigma)^2)
```

The program therefore sets:

```text
sigma = |y|
```

which produces an effective weighting proportional to:

```text
1 / y^2
```

The purpose of this weighting is to prevent the high-RFU standards from dominating the fit simply because their absolute residuals are numerically much larger than those of the low-RFU standards.

Because a zero or extremely small `sigma` would produce an excessively large weight, a small lower limit is applied to the sigma values.

The fitting function returns:

* Optimized 4PL parameters
* Parameter covariance matrix
* Estimated parameter errors
* Predicted RFUs
* Residuals
* Weighted residuals
* Sigma values used during fitting

---

# Standard-Curve Visualization

After fitting, the program generates a figure containing:

1. The observed standard measurements and fitted 4PL curve
2. The weighted residuals for the standards

The standard-curve plot uses logarithmic x and y axes. The residual plot uses a logarithmic x-axis.

The plot provides a visual method of evaluating whether:

* the standards follow the expected response,
* the fitted curve adequately describes the data, and
* individual standards show unusually large deviations from the model.

---

# Interpolation of Unknown Samples

Once the 4PL parameters have been determined, the equation is algebraically inverted to estimate nucleic-acid mass from an observed RFU.

The inverse relationship used by the program is:

```text
                              1 / Hill
             y - bottom
x = EC50 * (-------------)
               top - y
```

The implementation is:

```python
conc = ec50 * ((y-bottom)/(top-y)) ** (1/hill_slope)
```

This calculation converts:

```text
Measured RFU
    ->
Interpolated nucleic-acid mass in the quantification well
```

Values outside the fitted range are handled separately rather than blindly extrapolated.

---

# Dilution Correction

The interpolated value represents the amount of nucleic acid present in the Quantifluor measurement well, not the concentration of the original sample.

The program therefore applies process-specific dilution corrections.

Current dilution factors are:

| Sample Type | Dilution Correction |
| ----------- | ------------------: |
| FFPE DNA    |                  50 |
| FFPE RNA    |                  25 |
| cfDNA       |                   5 |
| gDNA        |                 100 |

The final concentration is calculated as:

```python
(interpolated_mass / 10) * dilution_factor
```

The division by 10 converts the interpolated mass to the concentration of the diluted sample because 10 µL of diluted sample was transferred into the Quantifluor assay.

The appropriate dilution correction is then applied to estimate the concentration of the original sample.

The resulting value is reported as:

```text
final_concentration (ng/ul)
```

---

# Reportability Checks

The program evaluates whether a calculated sample concentration can be reported.

Samples may be flagged when:

* the interpolated mass cannot be calculated,
* the fluorescence measurement is outside the fitted response range, or
* another invalid value prevents calculation of the final concentration.

Samples without an identified problem are assigned:

```text
reportable
```

This prevents invalid measurements from being treated as ordinary quantitative results.

---

# Batch-Level Quality Control

Several QC values are calculated for each quantification run.

## Positive Control Concentration

The calculated concentration of the designated control is extracted from the plate.

## Maximum Standard CV

The replicate %CV values are calculated for all standards, and the largest observed value is reported as:

```text
standard_max_cv
```

This provides a summary measure of the worst replicate agreement in the standard series.

## Maximum Standard Residual

The fitted standard curve is also used to back-calculate the mass of each standard.

The percent residual is calculated as:

```text
Interpolated Mass - Expected Mass
--------------------------------- x 100
          Expected Mass
```

The largest absolute residual is reported as:

```text
standard_max_residual
```

This provides a measure of how well the fitted calibration model reproduces the known standard concentrations.

---

# Process Validation

Before performing the final analysis, the program checks that the user's inputs agree with information contained in the Tecan file.

## Quantification Plate Barcode

The plate barcode entered by the user must match the plate barcode recorded in the uploaded `.asc` file.

If the two barcodes do not match, analysis stops.

This protects against accidentally analyzing a file from a different quantification plate.

## Quantifluor Method

The Tecan method embedded in the output file is mapped to an expected sample type:

```text
Quantifluor 2 DNA.mth   -> ffpe_dna
Quantifluor 2 RNA.mth   -> ffpe_rna
Quantifluor 2 cfDNA.mth -> cfdna
Quantifluor 2 gDNA.mth  -> gdna
```

The sample type selected by the user must agree with the method that generated the Tecan file.

This prevents the program from applying an incorrect standard configuration or dilution correction to the data.

---

# Output

The program generates two output formats.

## JSON

The JSON output is intended to provide structured data for downstream software.

It contains run-level information including:

```text
bottom_fit
top_fit
ec50_fit
hill_slope_fit
control_concentration
standard_max_cv
standard_max_residual
quant_data
```

The `quant_data` section contains the processed well-level quantification information.

## CSV

A CSV representation of the processed plate DataFrame is also generated.

The CSV is useful for:

* manual review,
* troubleshooting,
* development,
* opening results in Excel or similar software.

Both files are generated by the main workflow.

---

# Graphical User Interface

The program includes a Tkinter GUI to make the workflow accessible without requiring users to interact directly with Python.

The interface requests:

* Input `.asc` file
* Output file location
* Matrix plate barcode
* Quantification plate barcode
* Sample type

The available sample types are:

```text
ffpe_dna
ffpe_rna
cfdna
gdna
```

The user then selects:

```text
Run Quantification
```

The GUI validates required fields before running the analysis.

If the analysis succeeds, the status is changed to `Complete` and the user is informed where the results were saved.

If an exception occurs, the status is changed to `Failed` and the error message is displayed to the user.

---

# Major Program Components

The major functions in the program are:

### `query_wlo_extraction_plate()`

Queries Wet Lab Ops for the contents of a matrix plate.

### `_build_sigma_y2()`

Constructs the sigma values required to implement `1/y^2` weighted nonlinear regression.

### `four_param_logistic()`

Defines the four-parameter logistic calibration model.

### `four_pl_weighted_curve_fit()`

Fits the Quantifluor standard curve and calculates residuals and parameter uncertainty.

### `_plot_fit_and_residuals()`

Plots the fitted standard curve and weighted residuals.

### `parse_quant_file()`

Parses the Tecan `.asc` file, extracts metadata, maps plate positions, and retrieves tube barcodes from WLO.

### `four_pl_interpolate_concentration()`

Inverts the fitted 4PL equation to estimate nucleic-acid mass from fluorescence.

### `determine_concentration_status()`

Determines whether an individual sample concentration is reportable.

### `interpolate_raw_data()`

Performs the primary quantitative analysis, including standard processing, curve fitting, interpolation, dilution correction, and QC calculation.

### `run_quant_program()`

Coordinates the complete workflow and writes the output files.

### `submit()`

Receives input from the graphical interface, validates required fields, runs the analysis, and reports success or failure to the user.

---

# Summary

The Quantifluor 2 Quantification Program converts raw plate-reader fluorescence measurements into traceable nucleic-acid concentration results.

Rather than treating quantification as an isolated mathematical calculation, the program combines:

```text
Instrument output
       +
Sample identity
       +
Calibration modeling
       +
Dilution correction
       +
Quality control
       +
Process validation
       =
Structured quantification results
```

The purpose of the application is to reduce manual data manipulation, preserve sample identity, standardize quantification calculations, detect common process errors, and generate data in a format suitable for both human review and downstream automated processing.












example format:
ST1_1	A1	44307
ST1_2	B1	14857
ST1_3	C1	3646
ST1_4	D1	933
ST1_5	E1	265
ST1_6	F1	95
ST1_7	G1	54
BL1	H1	37
ST2_1	A2	45009
ST2_2	B2	14914
ST2_3	C2	3696
ST2_4	D2	965
ST2_5	E2	274
ST2_6	F2	100
ST2_7	G2	55
BL2	H2	40
SM1	A3	6559
SM4	B3	6030
SM7	C3	2670
SM10	D3	120
SM13	E3	515
SM16	F3	264
SM19	G3	8
SM22	H3	8
SM2	A4	7
SM5	B4	8
SM8	C4	9
SM11	D4	9
SM14	E4	9
SM17	F4	9
SM20	G4	8
SM23	H4	8
SM3	A5	8
SM6	B5	8
SM9	C5	8
SM12	D5	8
SM15	E5	8
SM18	F5	9
SM21	G5	8
SM24	H5	8
SM25	A6	8
SM28	B6	9
SM31	C6	8
SM34	D6	9
SM37	E6	8
SM40	F6	8
SM43	G6	8
SM46	H6	8
SM26	A7	8
SM29	B7	9
SM32	C7	9
SM35	D7	8
SM38	E7	8
SM41	F7	9
SM44	G7	10
SM47	H7	8
SM27	A8	8
SM30	B8	9
SM33	C8	8
SM36	D8	9
SM39	E8	8
SM42	F8	8
SM45	G8	8
SM48	H8	7
SM49	A9	9
SM52	B9	9
SM55	C9	8
SM58	D9	8
SM61	E9	9
SM64	F9	8
SM67	G9	8
SM70	H9	7
SM50	A10	8
SM53	B10	9
SM56	C10	8
SM59	D10	8
SM62	E10	8
SM65	F10	8
SM68	G10	7
SM71	H10	7
SM51	A11	8
SM54	B11	8
SM57	C11	8
SM60	D11	9
SM63	E11	9
SM66	F11	9
SM69	G11	8
PC1	H11	145
SM72	A12	8
SM73	B12	8
SM74	C12	8
SM75	D12	9
SM76	E12	8
SM77	F12	8
SM78	G12	8
PC2	H12	44	
Date of measurement: 2026-07-23/Time of measurement: 15:29:06
Quantifluor 2 DNA.mth
C:\Users\Public\Documents\Tecan\Magellan Pro\wsp\Buffy_Test1.wsp
490nm - 530nm
Scan Barcode: singlet_test
infinite 200Pro
Instrument serial number: 2106010740
Plate
Plate Description: [COS96fb] - Costar 96 Flat Black
Plate with Cover: No
Barcode: No
  Move Plate
  Action: Move out
  User Request
  Text: Feed the Tecan. Ca-caw.
  Move Plate
  Action: Move in
  Part of Plate
  Range: A1:H12
    Fluorescence Intensity
    Excitation Wavelength: 490 nm
    Excitation Bandwidth: 9 nm
    Emission Wavelength: 530 nm
    Emission Bandwidth: 20 nm
    ReadingMode: Top
    Lag Time: 0 �s
    Integration Time: 20 �s
    Number of Reads: 25
    Settle Time: 0 ms
    Gain: Optimal
    Z-Position: Manual
    Z-Position height: 20000 �m
    Label: test
Value of optimal gain: 77
Meas. temperature: Raw data: 20.5 �C
Date: 2026-07-23, Time: 15:29:06
