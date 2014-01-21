import xml.etree.cElementTree as ET

pepxml_file = "FIMC_digest_test04_1_100_pepsin-01.pep.xml"

for event, elem in ET.iterparse(pepxml_file):
    if elem.tag == '{http://regis-web.systemsbiology.net/pepXML}spectrum_query':
        print "Start"
        print elem[0][0].attrib['peptide']
        print elem.attrib['spectrum'].split(".")[0]
        print float(elem.attrib['retention_time_sec']) * 60
