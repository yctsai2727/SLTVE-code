function pdf=PDFnormalize(pdf,dx,dy)

temp=sum(pdf(:))*dx*dy;
pdf=pdf/temp;