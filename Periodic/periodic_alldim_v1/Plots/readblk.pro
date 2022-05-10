; (c) 2004 ANTI-RSI Foundation, and ANTI-IDL partnership.
; This sofware is given AS IS, and the author is not responsible 
; for any damage you can cause using this program.

pro readblk, file, n, out, NMB = nmb, NODES = nodes, ENDIAN = endian, $
        DIM = dim, ONLY = only, START = start, EXTN = extn

;+
; NAME:
;       READBLK
;
; PURPOSE:
;       Load binary files splitted in several blocks (as 
;       results of parallel runs) in 2D and 3D, in little 
;       endian or big endian format. 
;
; CATEGORY:
;       Input/output
;
; CALLING SEQUENCE:
;       READBLK, FILE, N, OUT
;
; INPUTS:
;       FILE:   The file name without numbers or extension
;
;       N:      Number of grid points
;
;       OUT:    Output variable
;
; KEYWORD PARAMETERS:
;       NMB:    The file number. If this keyword is not set, 
;               it will read *.in files
;
;       NODES:  The number of nodes used in the parallel run.
;
;       ENDIAN: If this keyword is set, the byte ordering of 
;               the array is reversed, making big endian into
;               little endian or vice-versa.
;
;       DIM:    The number of dimensions in the array
;
;       ONLY:   If this keyword is set, it will read only the 
;               outputs of the first 'ONLY' nodes.
;
;       START:  If ONLY and this keyword are set, it will read 
;               only the outputs of the first 'ONLY' nodes, 
;               starting at the node number 'START'. In any
;               other case this keyword is ignored.
;
;       EXTN:   File extension. If NODES is not set, the 
;               default extension is '.in'. In any other case
;               the extension is '.out' and this keyword is 
;               ignored.
;
; OUTPUTS:
;       The array is stored in OUT.
;
; COMMON BLOCKS:
;       None.
;
; RESTRICTIONS:
;       None. Anyone can use this stuff. However, you have to pay an 
;       IDL licence to RSI, to use the software that I made for free.
;       That is strange...
;
; AUTHOR:
;       Pablo Daniel Mininni (c)
;       Someone that wanted an IDL command 
;       to do the right thing ;-)
;-

on_error,2                      ;Return to caller if an error occurs
if keyword_set(nodes) eq 0 then nodes = 1
if keyword_set(dim) eq 0 then dim = 2
n1 = 1 & n2 = N & nprocs = nodes
if keyword_set(only) eq 0 then begin
   only = nodes
   start = 1
   nz = n
endif else begin
   count = 1
   if keyword_set(start) eq 0 then start = 1
   for i = 0,start+only-2 do begin
       irank = i
       iwork1 = (n2-n1+1)/nprocs
       iwork2 = n2-n1+1 MOD nprocs
       ista = irank*iwork1+n1+MIN(irank,iwork2)
       iend = ista+iwork1-1
       if (iwork2 gt irank) then begin 
          iend = iend+1
       endif
       ista = ista-count
       iend = iend-count
       if (i eq start-1) then sta = ista
       count = count+1
   endfor
   nz = iend-sta+1
endelse

if (dim eq 2) then OUT = fltarr(n,nz)
if (dim eq 3) then OUT = fltarr(n,n,nz)

trash = fltarr(1)

if keyword_set(nmb) eq 1 then begin
i = nmb                                   ;Generates file number string
if i le 9 then begin
   extn = strjoin([string(0,FORMAT='(I1)'), $
   string(0,FORMAT='(I1)'),string(i,FORMAT='(I1)'),'.out'])
endif
if (i ge 10) and (i le 99) then begin
   extn = strjoin([string(0,FORMAT='(I1)'), $
   string(i,FORMAT='(I2)'),'.out'])
endif
if (i ge 100) then begin
   extn = strjoin([string(i,FORMAT='(I3)'),'.out'])
endif
endif else begin
if keyword_set(extn) eq 0 then extn = 'in' ;Reads *.in files
endelse

count = 1                                 ;Reads files
for i = start-1,only-1 do begin
    irank = i
    iwork1 = (n2-n1+1)/nprocs
    iwork2 = n2-n1+1 MOD nprocs
    ista = irank*iwork1+n1+MIN(irank,iwork2)
    iend = ista+iwork1-1
    if (iwork2 gt irank) then begin 
       iend = iend+1
    endif
    ista = ista-count
    iend = iend-count
    count = count+1
    if i le 9 then begin
       node = strjoin(['.',string(0,FORMAT='(I1)'), $
       string(0,FORMAT='(I1)'),string(i,FORMAT='(I1)'),'.'])
    endif
    if (i ge 10) and (i le 99) then begin
       node = strjoin(['.',string(0,FORMAT='(I1)'), $
       string(i,FORMAT='(I2)'),'.'])
    endif
    if (i ge 100) then begin
       node = strjoin(['.',string(i,FORMAT='(I3)'),'.'])
    endif
    if (dim eq 2) then tmp = fltarr(n,iend-ista+1)
    if (dim eq 3) then tmp = fltarr(n,n,iend-ista+1)
    close,1
    openr,1,strjoin([FILE,node,extn])
    readu,1,trash
    readu,1,tmp
    close,1
    if (dim eq 2) then out(*,ista:iend) = tmp
    if (dim eq 3) then out(*,*,ista:iend) = tmp
endfor

if keyword_set(endian) eq 1 then swap_endian_inplace,out

end
