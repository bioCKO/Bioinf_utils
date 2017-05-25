__author__ = 'mjohnpayne'

import reportlab

from reportlab.pdfgen import canvas
from reportlab.lib.units import inch, cm
c = canvas.Canvas('/Volumes/MP_HD/CI_characterisation pics (possibly write display script?)/image_script/ex.pdf')
c.drawImage('/Volumes/MP_HD/CI_characterisation pics (possibly write display script?)/4D/c-b n-n/C-BSA_25-012-1.jpg', 100, 100, 10*cm, 10*cm)
c.showPage()
c.save()