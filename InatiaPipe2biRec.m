function [t,b] = InatiaPipe2biRec(dpipe,tpipe,hrec)
    %%%%%InatiaPipe2biRec(dpipe,tpipe,hrec)%%%%%%%%%%%%
    %直径d,肉厚tのパイプの断面二次モーメントに等しい、高さh、肉厚t、間隙bの２つの長方形の寸法を計算する
    %openVSPにはPipeの設定がないので、spar設定で頑張るときに使用する
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Ipipe = pi/64*((dpipe+tpipe/2)^4-(dpipe-tpipe/2)^4);
    t = Ipipe/2*12/hrec^3;
    b = sqrt((Ipipe/2-hrec*t^3/12)/hrec/t)*2;
end
