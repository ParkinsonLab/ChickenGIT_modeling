import os
import xml.etree.ElementTree as ET
import pandas as pd
import re
import numpy

# Function to extract EC numbers from an XML file
def extract_ec_numbers(xml_file):
    ec_numbers = set()
    tree = ET.parse(xml_file)
    root = tree.getroot()

    # Search for the pattern "ec-code:" followed by digits and dots in the XML file (in gapseq files)
    ec_pattern1 = re.compile(r"ec-code:(\d+\.\d+\.\d+\.\d+)")
    # or if that's cobrapy_adapted - then the pattern is "ec-code/X.X.X.X":
    ec_pattern2 = re.compile(r"ec-code/(\d+\.\d+\.\d+\.\d+)")

    # Combine matches from both patterns
    matches1 = ec_pattern1.findall(ET.tostring(root).decode())
    matches2 = ec_pattern2.findall(ET.tostring(root).decode())

    matches = matches1 + matches2
    for match in matches:
        ec_numbers.add(match)

    return ec_numbers

# Define the folder containing XML files
xml_folder_d10_models = "gapseq_models/2023/modelset_incl_gapseqAdapted_cobrapyAdapted_July2023/"

ec_numbers_models = {}

# Iterate over files with models (.xml files)
for file in os.listdir(xml_folder_d10_models):
    model_name = os.fsdecode(file)
    print(model_name)
    xml_file = os.path.join(xml_folder_d10_models,
                                model_name)
    ## THERE MIGHT BE AN ERROR THAT POPS UP - IF THERE ARE NON-XML FILEs (E.g. .DS.store)
    # make sure to read only .xml files in that case
    ec_numbers_in_xml = extract_ec_numbers(xml_file)
    # remove "_cobrapy_adapted.xml" and "-gapseqAdRm.xml" from model_name:
    model_name = model_name.replace("_cobrapy_adapted.xml", "").replace("-gapseqAdRm.xml", "")
    ec_numbers_models[model_name] = ec_numbers_in_xml

# Save the ec numbers for each model to a  file
# Open the file for writing
with open("ModelsAnalyses_outputs/Table_ECnumbers_inEachModel_gapseqAdapted_cobrapyAdapted_July2023.txt", 'w') as file:
    # Write the header
    file.write("SpecieName\tEC_numbers\n")

    # Iterate over each item in the dictionary
    for species, ec_numbers in ec_numbers_models.items():
        # Join the EC numbers into a comma-separated string
        ec_numbers_string = ', '.join(sorted(ec_numbers))
        # Write the species name and EC numbers to the file
        file.write(f"{species}\t{ec_numbers_string}\n")