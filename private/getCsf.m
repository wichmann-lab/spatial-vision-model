function [fgrid,sens] = getCsf(csfselector,timecourse)
%function [fgrid,sens] = getCsf(timecourse)
% a list of the csfs "fitted"
if strcmp(csfselector,'mean')
    if length(timecourse) <= 25
        fgrid = [0.05,.25,.5  ,1  ,2   ,3   ,4   ,8  ,16 ,32 ,64  ,100];
        sens  = [0.01,1.4,1.4 ,1.4,1.2 ,0.9 ,0.7 ,0.4,.3 ,.3 ,.001,0  ];
    elseif length(timecourse) <= 100
        fgrid = [0.05,.25 ,.5  ,1  ,2   ,3   ,4  ,8  ,12 ,16 ,32  ,64 ,100];
        sens  = [0.01,.5  ,1.3 ,2  ,1.9 ,1.9 ,1.8,1.6,1.3,1.2,.8  ,.1  ,0  ];
    elseif length(timecourse) <= 300
        fgrid = [0.05,.125,.25,.5 ,1 ,2 ,4  ,6  ,8 , 16,32  ,64  ,100];
        sens  = [0.01,0.01,.05, 3 ,5 ,5 ,6  ,7  ,8 , 8 ,2.5 ,.01 ,0  ];
    elseif length(timecourse) <=550
        fgrid = [0.05,.25,.5 ,1  ,2   ,4  ,6  ,8  , 16,32 ,64  ,100];
        sens  = [0.01,1  ,1.1,1.2,1.4 ,2.4,2.1,1.9,1 ,.7 ,.01 ,0  ];
    else
        fgrid = [0.05,.25,.5 ,1  ,2   ,4  ,8   ,16  ,32  ,64  ,100];
        sens  = [0.01, .3,.55,.8 ,1.4 ,2.2,1.65,1.25,1   ,.01 ,0  ];
    end
else
    if length(timecourse) <= 25
        fgrid = [0.05,.25,.5,1,2  ,3  ,4   ,8   ,16 ,32 ,64  ,100];
        sens  = [0.01,1.2,1.2 ,1.2,0.9,0.7 ,0.5 ,0.2,.1 ,.3 ,.001,0  ];
    elseif length(timecourse) <= 100
        fgrid = [0.05,.25 ,.5  ,1  ,2   ,3   ,4   ,8   ,11, 20 ,32  ,64 ,100];
        sens  = [0.01,0.05,0.4 ,.8 ,1.3 ,1.4 ,1.1 ,.8  ,.7, .6 ,.5  ,.2 ,0  ];
    elseif length(timecourse) <=550
        fgrid = [0.05,.25,.5 ,1  ,2   ,4  ,6  ,8  , 16,32 ,64  ,100];
        sens  = [0.01,1.5,1.7,2  ,2.4 ,2.4,2.1,1.9,1 ,.7 ,.01 ,0  ];
    else
        fgrid = [0.05,.25,.5 ,1  ,2   ,4  ,8 ,16  ,32  ,64  ,100];
        sens  = [0.01, .2,.4 ,.9 ,1.8 ,2.5,2 ,1.75,1   ,.01 ,0  ];
    end
end