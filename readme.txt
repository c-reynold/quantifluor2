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
