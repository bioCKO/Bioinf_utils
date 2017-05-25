__author__ = 'mjohnpayne'

text = open("/Users/mjohnpayne/PycharmProjects/Bioinf_utils/pythonchallenge/challenge3_text",'r').read().split('\n')
# print len(text)
# text = text.replace('\n',"")


length = "yYHToynVVpVjQFDHJmnJVjNuJIaoHwDvRGusLtSdWfmagxsaeJKHQHLCNCRiFTbpBLbofRNDpeLwEFOi"

print len(length)
import re



string = ''
for i in range(len(text)):
    if re.search("[A-Q]{3}[a-q][A-Q]{3}",text[i]):
        match = re.findall("[A-Q]{3}[a-q][A-Q]{3}",text[i])[0]
        pos = text[i].find(match)+3
        print "\n\n\n"
        print text[i-3][pos-4:pos+5]
        print text[i-2][pos-4:pos+5]
        print text[i-1][pos-4:pos+5]
        print text[i][pos-4:pos+5]
        print text[i+1][pos-4:pos+5]
        print text[i+2][pos-4:pos+5]
        print text[i+3][pos-4:pos+5]


# ans = re.findall("[a-q][A_Q]{2}[a-q][A-Q][a-q]",text)
#
#
#
# #ans = re.findall("[a-q][A-Q]{3}[a-q]",text)
# # out = ''
# # for i in ans:
# #     out += i[2:5]
# #
# # print out
# # #ans = re.findall("EXACTLY",text)
#
# print ans