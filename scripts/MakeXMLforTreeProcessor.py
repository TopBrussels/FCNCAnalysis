from shutil import copyfile
import xml.etree.cElementTree as ET
import os
cwd = os.getcwd()

#First copy the xml
copyfile(cwd+"/config/FullMcBkgdSamples_El.xml", cwd+"/config/FullMcBkgdSamples_El_TreeProcessor.xml")
copyfile(cwd+"/config/FullMcBkgdSamples_Mu.xml", cwd+"/config/FullMcBkgdSamples_Mu_TreeProcessor.xml")


channels = ["Mu","El"] 
# loop over channels
for chan in channels:
    print "\nSearching list of sample used for ", chan, " channel!"
    # getting the appropriate xml file
    if "Mu" in chan:
        filename = cwd+"/config/FullMcBkgdSamples_Mu_TreeProcessor.xml"
        tree = ET.ElementTree(file=cwd+'/config/FullMcBkgdSamples_Mu_TreeProcessor.xml')
    elif "El" in chan:
        filename = cwd+"/config/FullMcBkgdSamples_El_TreeProcessor.xml"
        tree = ET.ElementTree(file=cwd+'/config/FullMcBkgdSamples_El_TreeProcessor.xml')
    else:
        print "Channel '", chan , "' is not a correct channel name. No tree has been loaded!"
        sys.exit()

    root = tree.getroot()
    datasets = root.find('datasets')
    print "found  "  + str(len(datasets)) + " datasets"

    Luminosity_=0.;    

    for d in datasets:
        if "Data" in str(d.attrib['name']) and d.attrib['add'] == '1':
            Luminosity_=Luminosity_+float(d.attrib['EqLumi'])
            d.attrib['add'] = "0"

    # Append new tag: <a x='1' y='abc'>body text</a>
    new_tag = ET.SubElement(datasets, 'd')
    new_tag.attrib['name'] = "Data" # must be str; cannot be an int
    new_tag.attrib['title'] = "Data" # must be str; cannot be an int
    new_tag.attrib['add'] = "1"
    new_tag.attrib['color'] = "1"
    new_tag.attrib['ls'] = "1"
    new_tag.attrib['lw'] = "2"
    new_tag.attrib['normf'] = "1"
    new_tag.attrib['EqLumi'] = str(round(Luminosity_,9))
    new_tag.attrib['xsection'] =  "1."
    new_tag.attrib['PreselEff'] =  "0.0"
    new_tag.attrib['filenames'] =  "*.root"

    # Write back to file
    tree.write(filename)

#Next read all the lines
f_el = open(cwd+"/config/FullMcBkgdSamples_El_TreeProcessor.xml","r")
lines_el = f_el.readlines()
f_el.close()

#Remove the <data> and </data> lines
f_el = open(cwd+"/config/FullMcBkgdSamples_El_TreeProcessor.xml","w")
for line in lines_el:
  if line!="<data>"+"\n" and line!="</data>"+"\n":
    f_el.write(line)
f_el.close()
    
#Next read all the lines
f_mu = open(cwd+"/config/FullMcBkgdSamples_Mu_TreeProcessor.xml","r")
lines_mu = f_mu.readlines()
f_mu.close()

#Remove the <data> and </data> lines
f_mu = open(cwd+"/config/FullMcBkgdSamples_Mu_TreeProcessor.xml","w")
for line in lines_mu:
  if line!="<data>"+"\n" and line!="</data>"+"\n":
    f_mu.write(line)
f_mu.close()

