
__author__ = 'mjohnpayne'

text = open("/Users/mjohnpayne/PycharmProjects/Bioinf_utils/pythonchallenge/challenge2_text",'r').read()
print len(text)
text = text.replace('\n',"")

counts = {}
for i in text:
    if i not in counts:
        counts[i] = text.count(i)
        if counts[i] == 1:
            print i

#print counts