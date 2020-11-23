function baseline_beta = exp_baseline_fit(sfreqs,spect, verbose)

if nargin<3
    verbose=true;
end


%Set up data
y=prctile(spect(1:min(size(spect,1),20),:),1);
spect_bounds=prctile(y,[.1 99]);

w=sfreqs;

xl=[1.5 80];

% %Kill data bump and 60Hz noise
inds= (w<xl(1)) | (w>3 & w<35) | (w>58 & w<62)| (w>178 & w<182) |  (w>xl(2));

w=w(~inds);
y=y(~inds);


beta0=[abs(diff(spect_bounds)) .01 spect_bounds(1)];

model_fun=@(b,x)(b(1).*(1-b(2)).^x+b(3));

%Fit function
opts = statset('nlinfit');
opts.MaxIter=1000;
opts.TolX=1e-10; 

try
    baseline_beta = nlinfit(w,y,model_fun,beta0,opts);
catch
    warning('Did not properly fit baseline. Using default values');
    baseline_beta=beta0;
end

%Plot
if verbose
    figure
    plot(w,y,'linewidth',4)
    hold all;
    xhat=linspace(0,w(end),1000);
    plot(xhat,model_fun(baseline_beta,xhat),'linewidth',2);
    xlim(xl);
   
end
end
