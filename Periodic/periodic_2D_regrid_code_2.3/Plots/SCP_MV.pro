N=256
X=2*!pi*findgen(N)/N
FOR OUT_N=1,20 DO BEGIN
print,"FILE #",OUT_N
HEAD='./adv2Dps'
A_FN='a'
V_FN='v'
X_FN=['x','y','z']

FILE_AX=HEAD+A_FN+X_FN(0)
FILE_AY=HEAD+A_FN+X_FN(1)
FILE_AZ=HEAD+A_FN+X_FN(2)

print,FILE_AX
print,FILE_AY
print,FILE_AZ

print,"Reading A..."
readblk,HEAD,N,SS,dim=2,nmb=OUT_N,nodes=1

contour,SS,X,X,nlevels=32,/fill,/xs,/ys

ENDFOR

end
