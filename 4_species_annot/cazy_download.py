
import urllib2
import re
from bs4 import BeautifulSoup


import sys


records_per_page = 1000

prefix = "http://www.cazy.org/"

fivefamilies = ["Glycoside-Hydrolases.html","GlycosylTransferases.html","Polysaccharide-Lyases.html","Carbohydrate-Esterases.html","Carbohydrate-Binding-Modules.html", "Auxiliary-Activities.html"]
#fivefamilies = ["Carbohydrate-Binding-Modules.html"]

rx_member = re.compile(r'<a href="(\w+)\.html">')
rx_kingdom = re.compile(r'<a href="(http://www\.cazy\.org//\w+_(archaea|bacteria|eukaryota|viruses|characterized|structure)\.html)">\w+</a> (?:&#40;|\()(\d+).*(?:&#41;|\))</span>')
rx_subfamilies_exist = re.compile(r'http://www\.cazy\.org//\w+_subfamilies\.html">Subfamilies')
rx_ncbi = re.compile(r'http://www\.ncbi\.nlm\.nih\.gov/entrez/viewer\.fcgi\?db=protein\S+val=(\S+)"')
acc_cazy = {}

cazy_acc = {}

characterized = set()
structure = set()

done = set()
#go in to families
for family in fivefamilies:
    print family
    address = prefix + family


    f = urllib2.urlopen(address)

    html =  re.sub("</script />","</script>",f.read())

    html =  re.sub("<a href</td></tr>","<a href></td></tr>",html)
    html =  re.sub("<a href=G</td></tr>","<a href></td></tr>",html)

    soup = BeautifulSoup(html)

    for line in soup.findAll('td', width="30", align="right"):
        member = rx_member.search(str(line))
        if member:
            murl = prefix + member.group(1) + ".html"

            cazy_name = member.group(1)
            print cazy_name


            if murl not in done:
                done.add(murl)

                ff = urllib2.urlopen(murl)
                page = re.sub("</script />","</script>",ff.read())

                soup_member = BeautifulSoup(page)

                sub_family_exist = False
                if rx_subfamilies_exist.search(page):
                    sub_family_exist = True


                taxonurl = []
                for taxonline in soup_member.findAll('span'):

                    #print taxonline
                    taxon = rx_kingdom.search(str(taxonline))


                    #make a loop to collect all the urls first

                    if taxon:

                        taxonurl.append(taxon.group(1))


                        amount = int(taxon.group(3))
                        subpages_minus_one = (amount-1)/records_per_page

                        for i in xrange(subpages_minus_one):
                            taxonurl_address = prefix + "/" +  cazy_name + "_" + taxon.group(2) +  ".html?debut_PRINC=" + str((i+1)*1000) + "#pagination_PRINC"

                            taxonurl.append(taxonurl_address)

      #iterate the url here
                for taxonurl_address in taxonurl:

                    #print taxonurl_address
                    charac = False
                    structu = False


                    if  "characterized" in taxonurl_address:
                        charac = True
                    if  "structure" in taxonurl_address:
                        structu = True
                    taxonpage = urllib2.urlopen(taxonurl_address)
                    taxonpagecontent = re.sub("</script />","</script>",taxonpage.read())

                    taxon_member = BeautifulSoup(taxonpagecontent)


                    tabulka = taxon_member.findAll("tr", {"valign" : "top"})



                    for row in tabulka:


                        cols = row.findAll('td')


                        accession = ""
                        family_subfamily =cazy_name


                        search_ncbi = rx_ncbi.search(str(cols[3]))

                        if search_ncbi:
                            accession = search_ncbi.group(1).strip()



                        if sub_family_exist:

                            sub_family = str(BeautifulSoup(str(cols[6])).td.string)
                            if sub_family != "None":
                                family_subfamily +=   "-subfamily_" + str(BeautifulSoup(str(cols[6])).td.string)

                        if accession != "":
                            #print accession, family_subfamily

                            if accession not in acc_cazy:
                                acc_cazy[accession] = set()
                            acc_cazy[accession].add(family_subfamily)

                            if family_subfamily not in cazy_acc:
                                cazy_acc[family_subfamily] = set()
                            cazy_acc[family_subfamily].add(accession)

                            if charac == True:
                                characterized.add(accession)

                            if structu == True:
                                structure.add(accession)


#here we query ncbi and output the sequences

prefix = "http://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&rettype=fasta&id="
id_per_request = 200
path = sys.argv[1]

if not path.endswith("/"):
    path += "/"

for cazy_name in sorted(cazy_acc.keys()):
    print cazy_name
    temp_content = ""
    id_list = ""
    counter = 0
    temp_content = ""

    for acc in cazy_acc[cazy_name]:


        id_list += acc + ","
        counter += 1

        if counter == id_per_request:
            counter = 0
            url = prefix + id_list[:len(id_list)-1]
            #print url

            try:
                temp_content += urllib2.urlopen(url).read()

            except:
                for id in id_list[:len(id_list)-1].split(","):
                    url = prefix + id
                    #print url
                    try:
                        temp_content += urllib2.urlopen(url).read()
                    except:
                        print id
                        pass
            id_list = ""

    if id_list != "":
        url = prefix + id_list[:len(id_list)-1]
        id_list = ""

        try:
            temp_content += urllib2.urlopen(url).read()
        except:
            for id in id_list[:len(id_list)-1].split(","):
                url = prefix + id
                try:
                    temp_content += urllib2.urlopen(url).read()
                except:
                    print id
                    pass

    content = ""
    counts = 0
    for line in temp_content.splitlines():


        if ">" in line:
            counts += 1
            found = False
            content += line
            for acc in cazy_acc[cazy_name]:
                if acc in line.split("|"):

                    for cazy in acc_cazy[acc]:
                        content += "|" + cazy


                    if acc in characterized:
                        content += "_characterized"

                    if acc in structure:
                        content += "_structure"

                    found = True
                    break

            if found == False:
                print line + " no acc found in cazy_acc"
            content += "\n"
        else:
            content += line + "\n"
    print cazy_name + " expected:" +  str(len(cazy_acc[cazy_name])) + " real: " + str(counts)
    filename = path + cazy_name + ".fasta"
    f = open(filename, 'w')
    f.write(content)
    f.close()