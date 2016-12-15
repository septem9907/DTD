function [objval, gradient, hessian] = LikelihoodFunc(x, stage, data, h, maturity, initVal, adjVol, adjWgt)


%     if any(~isfinite(x))
%         objval = nan;
%         gradient = nan;
%         hessian = nan;
%         return
%     end

    switch stage
        case 1
            vol = x(1)*adjVol;
            wgt = x(2)*adjWgt;
        case 2
            vol = x(1)*adjVol;
            wgt = initVal.initWgt;
        otherwise
            % pass
    end
    
    % debt
    debt = data(:, h.curLiab)+0.5*data(:, h.ltBorrow)+wgt*data(:, h.othLiab);
    % implied asset
    impA = Solver(data, h, vol, debt, maturity, initVal.initImpA);
    n = length(impA);
    % Nd1
    sigsqrtt = vol*sqrt(maturity);
    d1 = (log(impA./debt)+(data(:, h.rfRate)+0.5*vol^2)*maturity)/sigsqrtt;
    
    % Nd1 = normcdf(d1);
    Nd1 = 0.5*erf(d1./sqrt(2)) + 0.5;
    Nd1(Nd1<eps) = eps;
    
    % objterm5
    numOfDate = sum(data(2:end, h.dateIntvl));
    logVA = log(impA./data(:, h.totAsset));
    logVAdif = diff(logVA);
    objterm5 = sum(logVAdif.^2./data(2:end, h.dateIntvl))-(logVA(end)-logVA(1))^2/numOfDate;
    % value of the objective function
    objval = (n-1)/2*log(2*pi)+(n-1)*log(vol)+sum(logVA(2:end))+sum(log(Nd1(2:end)))+0.5/vol^2*objterm5;
    
    
    if nargout > 1

                
        %% gradient        
        d2 = d1 - sigsqrtt;
                
        pd1 = normpdf(d1);
        pd2 = normpdf(d2);
        p2Nd1 = pd1./Nd1;
                
        exprt = exp(-data(:, h.rfRate)*maturity);
        expL2V = exprt.*debt./impA;
        expL2Vpd2 = expL2V.*pd2;
        
        mpd_Vsig_nu = pd1.*d2 - expL2Vpd2.*d1;
        mpd_Vsig_de = sigsqrtt*Nd1 + pd1 - expL2Vpd2;
        mpd_Vsig = sqrt(maturity) * (mpd_Vsig_nu ./ mpd_Vsig_de);   %main part of partial derivative : \frac{1}{V_t} \frac{\partial V_t}{\partial \sigma}
       
        mpd_Vsig_dif = mpd_Vsig(2:end)-mpd_Vsig(1:end-1);
        
        pd_d1sig = (mpd_Vsig/sqrt(maturity) - d2)/vol; % partial derivative: \frac{\partial d1}{\partial \sigma}
        
        grad1term3 = sum(p2Nd1(2:end).*pd_d1sig(2:end));
        grad1term5 = sum(logVAdif.*mpd_Vsig_dif./data(2:end, h.dateIntvl)) - (logVA(end)-logVA(1))*(mpd_Vsig(end)-mpd_Vsig(1))/numOfDate;
        
        % Calculate \frac{\partial L}{\partial \sigma}
        grad1 = (n-1)/vol + sum(mpd_Vsig(2:end)) + grad1term3 - objterm5/vol^3 + grad1term5/vol^2;
        
        if stage == 1
            c2L = data(:, h.othLiab)./debt;
            % Nd2 = normcdf(d2);
            Nd2 = 0.5*erf(d2./sqrt(2)) + 0.5;
            Nd2(Nd2<eps) = eps;
            
            mpd_Vdel_nu = pd1 - expL2Vpd2 + sigsqrtt*expL2V.*Nd2;
            mpd_Vdel_de = mpd_Vsig_de;
            mpd_Vdel = c2L .* (mpd_Vdel_nu ./ mpd_Vdel_de);   %main part of partial derivative : \frac{1}{V_t} \frac{\partial V_t}{\partial \delta}
            mpd_Vdel_dif = mpd_Vdel(2:end)-mpd_Vdel(1:end-1);
            pd_d1del = (mpd_Vdel - c2L)/sigsqrtt;  % partial derivative: \frac{\partial d1}{\partial \delta}
            
            grad2term2 = sum(p2Nd1(2:end).*pd_d1del(2:end));
            grad2term3 = sum(logVAdif.*mpd_Vdel_dif./data(2:end, h.dateIntvl)) - (logVA(end)-logVA(1))*(mpd_Vdel(end)-mpd_Vdel(1))/numOfDate;
            
            % Calculate \frac{\partial L}{\partial delta}
            grad2 = sum(mpd_Vdel(2:end)) + grad2term2 + grad2term3/vol^2;
            
            gradient = [grad1*adjVol; grad2*adjWgt];
        else
            gradient = grad1*adjVol;
        end
                
        
        %% hessian                       
        pd_d2sig = (mpd_Vsig/sqrt(maturity) - d1)/vol; % partial derivative: \frac{\partial d2}{\partial \sigma}
        
        
        msecpd_Vsigsig_1 = pd1.*(-d1.*d2.*pd_d1sig + pd_d2sig) - expL2Vpd2.*(-mpd_Vsig.*d1 - d1.*d2.*pd_d2sig + pd_d1sig);
        msecpd_Vsigsig_2 = sqrt(maturity)*Nd1 + pd1.*pd_d1sig.*(sigsqrtt-d1) - expL2Vpd2.*(-mpd_Vsig - d2.*pd_d2sig);
        
        % Calculate -\frac{1}{V_t^2}(\frac{\partial V_t}{\partial \sigma})^2 + \frac{1}{V_t} \frac{\partial^2 V_t}{\partial \sigma^2}
        msecpd_Vsigsig =  sqrt(maturity)./mpd_Vsig_de.^2 .* (msecpd_Vsigsig_1.* mpd_Vsig_de - mpd_Vsig_nu.* msecpd_Vsigsig_2);
                
        secpd_d1sigsig = (-pd_d1sig + (msecpd_Vsigsig/sqrt(maturity) - pd_d2sig))/vol; % second partial derivative: \frac{\partial^2 d1}{\partial \sigma^2}
        
        msecpd_Vsigsig_dif = msecpd_Vsigsig(2:end) - msecpd_Vsigsig(1:end-1);

        hessian11term3vec = p2Nd1.*(-d1.*pd_d1sig.^2 + secpd_d1sigsig) - p2Nd1.^2.*pd_d1sig.^2;
        hessian11term6 =  sum((mpd_Vsig_dif.^2 + logVAdif.* msecpd_Vsigsig_dif)./data(2:end, h.dateIntvl))...
            - ((mpd_Vsig(end)-mpd_Vsig(1))^2 + (logVA(end)-logVA(1))*(msecpd_Vsigsig(end)-msecpd_Vsigsig(1)))/numOfDate;
        
        % Calculate \frac{\partial^2 L}{\partial \sigma^2}                
        hessian11 = -(n-1)/vol^2 + sum(msecpd_Vsigsig(2:end)) + sum(hessian11term3vec(2:end)) + 3/vol^4*objterm5 - 4/vol^3*grad1term5 + hessian11term6/vol^2;
                
        if stage == 1
            pd_d2del = pd_d1del;  % partial derivative: \frac{\partial d2}{\partial \delta}
            expmix = exprt.*(data(:, h.othLiab)./impA - debt./impA.*mpd_Vdel);
            pNd2dif = pd2 - sigsqrtt*Nd2;
            
            msecpd_Vdeldel_1 = -pd1.*d1.*pd_d1del -  expmix.*pNd2dif - expL2Vpd2.*pd_d2del.*(-d2-sigsqrtt);
            msecpd_Vdeldel_2 = pd1.*pd_d1del.*(sigsqrtt-d1) - expmix.*pd2 + expL2Vpd2.*d2.*pd_d2del;
            
            % Calculate -\frac{1}{V_t^2}(\frac{\partial V_t}{\partial \delta})^2 + \frac{1}{V_t} \frac{\partial^2 V_t}{\partial \delta^2}
            msecpd_Vdeldel = c2L.*(-mpd_Vdel + (msecpd_Vdeldel_1.*mpd_Vdel_de - mpd_Vdel_nu.*msecpd_Vdeldel_2)./mpd_Vdel_de.^2);
            
            msecpd_Vsigdel_1 = -pd1.*d1.*pd_d1sig + expL2V.* mpd_Vsig.*pNd2dif - expL2V.*(pd2.*pd_d2sig.*(-d2-sigsqrtt) - Nd2*sqrt(maturity));
            msecpd_Vsigdel_2 = msecpd_Vsigsig_2;
            
            % Calculate -\frac{1}{V_t^2}(\frac{\partial V_t}{\partial \sigma})*(\frac{\partial V_t}{\partial \sigma}) + \frac{1}{V_t} \frac{\partial^2 V_t}{\partial \sigma \partial \delta}
            msecpd_Vsigdel = c2L./mpd_Vdel_de.^2 .* (msecpd_Vsigdel_1.*mpd_Vdel_de - mpd_Vdel_nu.*msecpd_Vsigdel_2);
            
            % expmix2 = exprt.*(data(:,4)./evalue - L2V.* mpd_Vdel);
            %
            % msecpd_Vdelsig_1 = pd1.*(-d1.*d2.*pd_d1del + pd_d2del) - expmix2.*pd2.*d1 - expL2Vpd2.*(-d1.*d2.*pd_d2del + pd_d1del);
            % msecpd_Vdelsig_2 = msecpd_Vdeldel_2;
            %
            % Calculate -\frac{1}{V_t^2}(\frac{\partial V_t}{\partial \sigma})*(\frac{\partial V_t}{\partial \sigma}) + \frac{1}{V_t} \frac{\partial^2 V_t}{\partial \delta \partial \sigma}
            % msecpd_Vdelsig = sqrt(tdays)./mpd_Vsig_de.^2.*(msecpd_Vdelsig_1.*mpd_Vsig_de - mpd_Vsig_nu.*msecpd_Vdelsig_2);
            
            secpd_d1deldel = (msecpd_Vdeldel + c2L.^2)/sigsqrtt; % second partial derivative: \frac{\partial^2 d1}{\partial \delta^2}
            secpd_d1sigdel = -pd_d1del/vol + msecpd_Vsigdel/sigsqrtt;  % second partial derivative: \frac{\partial^2 d1}{\partial \sigma \partial \delta}
            % secpd_d1delsig = (msecpd_Vsigdel/sqrt(tdays) - pd_d2del)/sig; % second partial derivative: \frac{\partial^2 d1}{\partial \delta \partial \sigma}
            
            msecpd_Vdeldel_dif = msecpd_Vdeldel(2:end) - msecpd_Vdeldel(1:end-1);
            msecpd_Vsigdel_dif = msecpd_Vsigdel(2:end) - msecpd_Vsigdel(1:end-1);
            
            hessian22term2vec = p2Nd1.*(-d1.*pd_d1del.^2 + secpd_d1deldel) - p2Nd1.^2.*pd_d1del.^2;
            hessian22term3 =  sum((mpd_Vdel_dif.^2 + logVAdif.* msecpd_Vdeldel_dif)./data(2:end, h.dateIntvl))...
                - ((mpd_Vdel(end)-mpd_Vdel(1))^2 + (logVA(end)-logVA(1))*(msecpd_Vdeldel(end)-msecpd_Vdeldel(1)))/numOfDate;
            
            % Calculate \frac{\partial^2 L}{\partial \delta^2}
            hessian22 = sum(msecpd_Vdeldel(2:end)) + sum(hessian22term2vec(2:end)) + hessian22term3/vol^2;
            
            hessian12term2vec = p2Nd1.*(-d1.*pd_d1sig.*pd_d1del + secpd_d1sigdel) - p2Nd1.^2.*pd_d1sig.*pd_d1del;
            hessian12term4 =  sum((mpd_Vsig_dif.*mpd_Vdel_dif + logVAdif.* msecpd_Vsigdel_dif)./data(2:end, h.dateIntvl))...
                - ((mpd_Vsig(end)-mpd_Vsig(1))*(mpd_Vdel(end)-mpd_Vdel(1)) + (logVA(end)-logVA(1))*(msecpd_Vsigdel(end)-msecpd_Vsigdel(1)))/numOfDate;
            
            % Calculate \frac{\partial^2 L}{\partial \sigma \partial \delta}
            hessian12 = sum(msecpd_Vsigdel(2:end)) + sum(hessian12term2vec(2:end)) - 2/vol^3*grad2term3 + hessian12term4/vol^2;
            
            hessian = [hessian11*adjVol^2, hessian12*adjVol; hessian12*adjVol, hessian22];            
        else
            hessian = hessian11*adjVol^2;
        end
        
        
    end

end