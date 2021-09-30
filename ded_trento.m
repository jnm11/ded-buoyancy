

nm='gc/f6/063';

%rsync -vaP --append-verify tpos:gc/f6/063 $DEDALUS_DATA/gc/pos/f6 --exclude "b*"
%rsync -vaP --append-verify asahi:gc/pos/f6/063 $DEDALUS_DATA/gc/pos/f6 --exclude "b*"
