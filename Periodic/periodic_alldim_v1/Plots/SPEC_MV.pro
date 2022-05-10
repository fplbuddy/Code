LOADCT,6
b0=5

NN=512/2
MM=512/2
X=findgen(MM/2+1)+1
Y=X

UU=dblarr(MM/2+1,NN/2+1)
BB=UU
KK=UU
for i=0,MM/2 DO BEGIN
for j=0,MM/2 DO BEGIN
KK(i,j)=( (i+1.0d)^2+(j+1.0d)^2 )
endfor
endfor

DIR='./'

FOR IT =1,100,1 do begin
                   ID=      strtrim(string(IT),1)
if (IT le 99) then ID='0'  +strtrim(string(IT),1) 
if (IT le 9 ) then ID='00' +strtrim(string(IT),1)

file1=DIR+'spectrum_UU.'+ID+'.out'
print,file1
close,1 & openr,1,file1
readu,1,tmp
readu,1,UU
close,1 
file1=DIR+'spectrum_BB.'+ID+'.out'
print,file1
close,1 & openr,1,file1
readu,1,tmp
readu,1,BB
close,1

;contour,ALOG(WW+1.0e-14)   ,X,Y,/xlog,/ylog,/fill,/xs,/ys,nlevels=32
contour,ALOG((UU+BB)*KK+1.0e-14)   ,X,Y,/xlog,/ylog,/fill,/xs,/ys,nlevels=32
oplot,[1,5,5],[5,5,1]
wait,0.5

ENDFOR

end
