./finalAnalysis.sh output/$1PbPb0_30 1_2  $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} --MinZPT 40 --MaxZPT 350 -MinTrackPT 1 --MaxTrackPT 2  --MinHiBin 0 --MaxHiBin 60  --Input mergedSample/HISingleMuon-v6.root --MixFile mergedSample/HISingleMuon-v6.root  --IsPP false --IsGenZ false --nMix 10
./finalAnalysis.sh output/$1PbPb0_30 2_4  $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} --MinZPT 40 --MaxZPT 350 -MinTrackPT 2 --MaxTrackPT 4  --MinHiBin 0 --MaxHiBin 60  --Input mergedSample/HISingleMuon-v6.root --MixFile mergedSample/HISingleMuon-v6.root  --IsPP false --IsGenZ false --nMix 10
./finalAnalysis.sh output/$1PbPb0_30 4_10 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} --MinZPT 40 --MaxZPT 350 -MinTrackPT 4 --MaxTrackPT 10 --MinHiBin 0 --MaxHiBin 60  --Input mergedSample/HISingleMuon-v6.root --MixFile mergedSample/HISingleMuon-v6.root  --IsPP false --IsGenZ false --nMix 10

hadd -f plots/$1PbPb0_30-result.root output/$1PbPb0_30-1_2-result.root output/$1PbPb0_30-2_4-result.root output/$1PbPb0_30-4_10-result.root

./finalAnalysis.sh output/$1PbPb30_50 1_2  $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} --MinZPT 40 --MaxZPT 350   -MinTrackPT 1 --MaxTrackPT 2  --MinHiBin 60 --MaxHiBin 100  --Input mergedSample/HISingleMuon-v6.root --MixFile mergedSample/HISingleMuon-v6.root  --IsPP false --IsGenZ false --nMix 10
./finalAnalysis.sh output/$1PbPb30_50 2_4  $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} --MinZPT 40 --MaxZPT 350  --MinTrackPT 2 --MaxTrackPT 4  --MinHiBin 60 --MaxHiBin 100  --Input mergedSample/HISingleMuon-v6.root --MixFile mergedSample/HISingleMuon-v6.root  --IsPP false --IsGenZ false --nMix 10
./finalAnalysis.sh output/$1PbPb30_50 4_10 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} --MinZPT 40 --MaxZPT 350  --MinTrackPT 4 --MaxTrackPT 10 --MinHiBin 60 --MaxHiBin 100  --Input mergedSample/HISingleMuon-v6.root --MixFile mergedSample/HISingleMuon-v6.root  --IsPP false --IsGenZ false --nMix 10

hadd -f plots/$1PbPb30_50-result.root output/$1PbPb30_50-1_2-result.root output/$1PbPb30_50-2_4-result.root output/$1PbPb30_50-4_10-result.root

./finalAnalysis.sh output/$1PbPb50_90 1_2  $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11}  --MinZPT 40 --MaxZPT 350 -MinTrackPT 1 --MaxTrackPT 2  --MinHiBin 100 --MaxHiBin 180  --Input mergedSample/HISingleMuon-v6.root --MixFile mergedSample/HISingleMuon-v6.root  --IsPP false --IsGenZ false --nMix 10
./finalAnalysis.sh output/$1PbPb50_90 2_4  $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11}  --MinZPT 40 --MaxZPT 350 -MinTrackPT 2 --MaxTrackPT 4  --MinHiBin 100 --MaxHiBin 180  --Input mergedSample/HISingleMuon-v6.root --MixFile mergedSample/HISingleMuon-v6.root  --IsPP false --IsGenZ false --nMix 10
./finalAnalysis.sh output/$1PbPb50_90 4_10 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11}  --MinZPT 40 --MaxZPT 350 -MinTrackPT 4 --MaxTrackPT 10 --MinHiBin 100 --MaxHiBin 180  --Input mergedSample/HISingleMuon-v6.root --MixFile mergedSample/HISingleMuon-v6.root  --IsPP false --IsGenZ false --nMix 10

hadd -f plots/$1PbPb50_90-result.root output/$1PbPb50_90-1_2-result.root output/$1PbPb50_90-2_4-result.root output/$1PbPb50_90-4_10-result.root

./finalAnalysis.sh output/$1PbPb0_90 1_2   $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} --MinZPT 40 --MaxZPT 350 -MinTrackPT 1 --MaxTrackPT 2  --MinHiBin 0 --MaxHiBin 180  --Input mergedSample/HISingleMuon-v6.root --MixFile mergedSample/HISingleMuon-v6.root  --IsPP false --IsGenZ false --nMix 10
./finalAnalysis.sh output/$1PbPb0_90 2_4   $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} --MinZPT 40 --MaxZPT 350 -MinTrackPT 2 --MaxTrackPT 4  --MinHiBin 0 --MaxHiBin 180  --Input mergedSample/HISingleMuon-v6.root --MixFile mergedSample/HISingleMuon-v6.root  --IsPP false --IsGenZ false --nMix 10
./finalAnalysis.sh output/$1PbPb0_90 4_10  $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} --MinZPT 40 --MaxZPT 350 -MinTrackPT 4 --MaxTrackPT 10 --MinHiBin 0 --MaxHiBin 180  --Input mergedSample/HISingleMuon-v6.root --MixFile mergedSample/HISingleMuon-v6.root  --IsPP false --IsGenZ false --nMix 10

hadd -f plots/$1PbPb0_90-result.root output/$1PbPb0_90-1_2-result.root output/$1PbPb0_90-2_4-result.root output/$1PbPb0_90-4_10-result.root
