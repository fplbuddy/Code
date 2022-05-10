LOADCT,6
b0=5

NN=512/2
MM=512/2
X=findgen(MM)/float(MM)*2*!pi*2
Y=X

WW=dblarr(MM,NN)
AA=WW
DIR='./'

FOR IT =1,100,1 do begin
                   ID=      strtrim(string(IT),1)
if (IT le 99) then ID='0'  +strtrim(string(IT),1) 
if (IT le 9 ) then ID='00' +strtrim(string(IT),1)

file1=DIR+'mhd2Dww.000.'+ID+'.out'
close,1 & openr,1,file1
readu,1,tmp
readu,1,WW
close,1 
file2=DIR+'mhd2Daa.000.'+ID+'.out'
close,1 & openr,1,file2
readu,1,tmp
readu,1,AA
close,1
print,file2

for i=0,NN-1 DO BEGIN
for j=0,NN-1 DO BEGIN
AA(i,j)=AA(i,j)+b0*(j-NN/2)*2*!pi/NN
endfor
endfor


contour,WW   ,X,Y,/fill,/xs,/ys,nlevels=32
contour,AA   ,X,Y,/overplot,levels=[-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]*max(AA)
;contour,AA   ,X,Y,/overplot,levels=[-0.4,0.0,0.4]*max(AA),thick=4
;print,min(WW),max(WW)
wait,0.5


ENDFOR

end
