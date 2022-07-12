function [Final_indir_EF_sec_per_hh, FinIND_EF_sec_hh_12COICOP, Final_indir_EF_hh, Final_EF_sec_per_hh, Final_Direct_EF_sec_hh, Final_Direct_EF_hh, Final_EF_sec_per_hh_12COICOP, Final_EF_sec_per_hh_47COICOP, Final_exp_sec_hh_tot, Final_EXP_sec_per_hh_12COICOP, Final_EXP_sec_per_hh_47COICOP] = Vietnam_EF_function()

%% Calculating direct and indirect energy for Nepal based on expenditures   

% 140 countries
ctry = 19; % choose here for which country you do your analysis
% 57 sectors
% Nepal is the 23th country in IEA and GTAP
% 23 IEA sectors

%%
%%load data
load GTAP_basics.mat;
IEA_to_GTAP_mask(:,:,1) = xlsread('O:/eemba/Data/IEA/GTAP_IEA_final_energy_extension_Yannick.xlsx','IEA_to_GTAP_dir_and_indir','C3:BG25');
IEA_FC_tot(:,:,1) = xlsread('O:/eemba/Data/IEA/IEA_world clean 2011_Yannick.xlsx','MB','P2:P3221'); %including bunkers and road indirect
IEA_road_dir (:,:,1) = xlsread('O:/eemba/Data/IEA/IEA_world clean 2011_Yannick.xlsx','MB','U2:U141');
IEA_resid_dir (:,:,1) = xlsread('O:/eemba/Data/IEA/IEA_world clean 2011_Yannick.xlsx','MB', 'T2:T141');


   
 %% creating the energy vector from GTAP industries shares and IEA to GTAP mapping
 
 %STEP I - loaded IEA to GTAP mask - already done before (look above:
 %IEA_to_GTAP_mask - it is matrix full of ones - many to many. 
           
 %STEPII: creating total spends on products per each industry for each country per sector (140x57)
 %take spends on on energy products by each industry. I choose
%following energy products from GTAP: 15 - coal, 16-oil, 17-gas,43-
%electricity, 44 - gas manufacture. I use here Z matrix
 Z = A .* repmat(X,1,7980);
sum_prod_per_ind  = sum(Z(15:57:end,:))+sum(Z(16:57:end,:))+sum(Z(17:57:end,:))+sum(Z(43:57:end,:))+sum(Z(44:57:end,:));

%STEP III: energy extension based on shares of spends

for c = 1:140%for countries from 1 to 140
     mask_times_spends_per_sec_within_country (:,(c*57)-56:(c*57)) = IEA_to_GTAP_mask .* repmat ((sum_prod_per_ind (:,(c*57)-56:(c*57))),23,1);
     spends_shares_sum_rowwise (:,(c*57)-56:(c*57))= mask_times_spends_per_sec_within_country (:,(c*57)-56:(c*57)) ./ repmat ((sum(( mask_times_spends_per_sec_within_country (:,(c*57)-56:(c*57))),2)),1,57);
     IEA_FC_tot_trans = IEA_FC_tot';
     Big_tot_IEA_to_GTAP (:,(c*57)-56:(c*57)) = IEA_FC_tot_trans (:,(c*23)-22:(c*23)) * spends_shares_sum_rowwise (:,(c*57)-56:(c*57));
 end
   
 %STEP IV: creating direct energy vectors: transport and residential:
%already loaded - look above for IEA_resid_dir and IEA_road_dir. Both are
%140x1 vectors. 

%% RAS BALANCING

% For RAS balancing you need to first have sum for GTAP Y_h for Nepal and
% sum of expenditures from Nepalese CES 

 % STEP 1: Take totals for GTAP-for example for agriculture (first sector) it will
 % be total of agr in australia, germany, uk etc. And totals for COICOP
 % which is a sum along the rows.
 
 for s =1:57
  GTAP_Viet_sec_tot (:,s) = sum(Y_h (s:57:end,ctry));
 end
 
 
 %% 
 % IMPORTANT! Depending on country you use for your analysis you should
 % input here specific data (load) and specify it within variable ctry_exp
 % load all of the expenditures for all households

 %HOUSEHOLDS WITHOUT 9 STD DEV OUTLIERS
Viet_exp_w_weights =  xlsread('O:/eemba/Data/Vietnam/yearlyEXP_VIET.xls','WGHT_Viet_yr_exp_wwout_9std_out','D2:EG8838');

Viet_exp_wghts_doll = ((Viet_exp_w_weights .*0.000051198).* 1.032  ) ./1000;

%multiplied with 1.032 to account for inflation - the same as for Nepal since it was the same period,
%and 0.000051198 to swich from
%vietnamese dong to US dollar (exchange rate for 2010-2011) and divided by thousand to get to
%the same format as GTAP has. In the Vietnamese survey the unit is thousand vietnamese dong!


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
Ctry_exp = Viet_exp_wghts_doll ; %choose data
numb_coicop_cat = size (Ctry_exp,2); % indicate how many expenditure categories your data has
numb_hh = size (Ctry_exp,1); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Continuation of step 1 for RAS Balancing
 %TOTAL for expenditures in Ctry_exp. 
 
Ctry_exp_total = (sum (Ctry_exp,1))';%140 x 1 TRANSPOSED!

 %check the difference between total for Ctry_exp_total and GTAP_Nep_sec_tot
 diff_GTAP_COICOP = sum (GTAP_Viet_sec_tot,2) - sum (Ctry_exp_total,1); 
 
 %make shares for Ctry_exp_total- what share each sector has in regards to total
 Ctry_exp_total_shares = Ctry_exp_total ./ repmat((sum(Ctry_exp_total,1)),numb_coicop_cat,1);
  %make a rebalanced total COICOP so the total for COICOP and GTAP_Nep_sec_tot is the same  
 Ctry_exp_total_rebal = Ctry_exp_total + (Ctry_exp_total_shares .* repmat (diff_GTAP_COICOP,numb_coicop_cat,1));

 
 %STEP 2: using mask COICOP_to_GTAP (zeros and ones) rebalance the matrix using RAS
%method to populate the inside of the matrix -instead of zeros and ones use
%numbers. So the sum of each column will be equal to the sum of rows. 

%load COICOP_to_GTAP masks:

VietCOICOPtoGTAP_mask =xlsread('O:/eemba/Data/Vietnam/Viet_COICOP_to_GTAP.xlsx','COICOP_to_GTAP_2010','F3:BJ136');

%first make sures you don't have any zeros in your totals
Ctry_exp_total_rebal(Ctry_exp_total_rebal<=0) = 0.00000001;
GTAP_Viet_sec_tot (GTAP_Viet_sec_tot <=0) = 0.00000001;
 n = 1000; %number of iterations
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
 RASbal = VietCOICOPtoGTAP_mask; % here I store final rebalanced matrix after RAS
 %and also depending on a number of expenditures in each country you will use different
 %COICO to GTAP mask
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%

for i = 1:n
    
  rowsumnow = sum (RASbal,2);   
  rowsumnow(rowsumnow<=0)=0.000000001; 
  
     for r= 1:numb_coicop_cat %for all rows corresponding to Vietnamese exp categories (134)
         RASbal (r,:) = RASbal (r,:) ./ rowsumnow (r) .* Ctry_exp_total_rebal (r);  
      end
  
  colsumnow = sum (RASbal,1);
  colsumnow(colsumnow<=0)=0.000000001;
    
     for c=1:57 %for 57 GTAP sectors
       RASbal (:,c) = RASbal (:,c) ./ colsumnow (:,c) .* GTAP_Viet_sec_tot (:,c);   
     end  
    
end

%% Stretch final demand to get to 140 coicop categoires for each household
%Y_h(Y_h<=0) = 0.00000001; %not used here
Y_h_Viet_shares = Y_h (:,ctry) ./ (repmat (GTAP_Viet_sec_tot',140,1));%if you sum every GTAP sector 1 for each country you will get 100%
Y_all_sec = repmat (RASbal',140,1) .* repmat (Y_h_Viet_shares,1,numb_coicop_cat); %7980x134
%% ENERGY FOOTPRINT
X(X==0) = 0.000000001; % the total output cannot have zeros
%A = Z./X';
I = eye(7980);
% L  = inv (I -A); L is already calculated
Xdiag = diag (X);
e = Big_tot_IEA_to_GTAP * inv (Xdiag) ; %energy intensities vector
Ctry_exp_shares_100perc_for_sec = Ctry_exp'./ repmat((sum(Ctry_exp',2)),1,numb_hh); %134x9399 matrix: the sum along rows (,2) gives 100%
 for cat = 1: numb_coicop_cat
   spends_per_sec_tmp = Y_all_sec (:,cat) * Ctry_exp_shares_100perc_for_sec (cat,:); %7980x9399 tot spends of all hh on the expenditure category cat
   Final_indir_EF_sec_per_hh(cat,:) = e*L *spends_per_sec_tmp ; %134x9399 matrix
 end
 
Final_indir_EF_hh=(sum(Final_indir_EF_sec_per_hh,1))'; 

% calulating footprint of each expenditure category for all hh
% Footprints for sectors related to direct energy use need to be calculated
% seperatly and added to this calculation 
 
 %%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %DIRECT ENERGY FOOTPRINTS
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%
 %Residential direct EF

%HOUSEHOLDS WITHOUT 9 STD DEV OUTLIERS
resid_exp=(((xlsread('O:/eemba/Data/Vietnam/yearlyEXP_VIET.xls','WGHT_Direct_exp_Viet_wwout_9std','D3:K8839')).*0.000051198).* 1.032) ./1000;
tot_resid_exp=sum(resid_exp,1);

IEA_resid_split = [(0.07*IEA_resid_dir(ctry,:)),(0.07*IEA_resid_dir(ctry,:)),(0.054*IEA_resid_dir(ctry,:)),(0.054*IEA_resid_dir(ctry,:)),(0.054*IEA_resid_dir(ctry,:)),(0.695*IEA_resid_dir(ctry,:)),(0.695*IEA_resid_dir(ctry,:)),(0.181*IEA_resid_dir(ctry,:))];

sum_resid_sec_coal=tot_resid_exp(:,1)+tot_resid_exp(:,2) ;
sum_resid_sec_oil = tot_resid_exp(:,3)+tot_resid_exp(:,4)+tot_resid_exp(:,5) ;
sum_resid_sec_biofuel =tot_resid_exp(:,6)+tot_resid_exp(:,7);

shares_resid_sec_coal=[(tot_resid_exp(:,1)/sum_resid_sec_coal), (tot_resid_exp(:,2)/sum_resid_sec_coal)];
shares_resid_sec_oil = [(tot_resid_exp(:,3)/sum_resid_sec_oil), (tot_resid_exp(:,4)/sum_resid_sec_oil), (tot_resid_exp(:,5)/sum_resid_sec_oil)];
shares_resid_sec_biofuel = [(tot_resid_exp(:,6)/sum_resid_sec_biofuel), (tot_resid_exp(:,7)/sum_resid_sec_biofuel)];

tot_IEA_resid=IEA_resid_split.*[shares_resid_sec_coal,shares_resid_sec_oil, shares_resid_sec_biofuel,1];

for r=1:8
shares_resid_exp(:,r)=resid_exp(:,r)./ repmat((sum(resid_exp(:,r),1)),numb_hh,1);
EF_resid_dir(:,r)= shares_resid_exp(:,r).*repmat((tot_IEA_resid(:,r)),numb_hh,1);
end

 %% calculating for direct road
 
%HOUSEHOLDS WITHOUT 9 STD DEV OUTLIERS   
dir_road_exp=(((xlsread('O:/eemba/Data/Vietnam/yearlyEXP_VIET.xls','WGHT_Direct_exp_Viet_wwout_9std','M3:O8839')).*0.000051198).* 1.032) ./1000;
tot_dir_road_exp=sum(dir_road_exp,1);
shares_dir_road=[tot_dir_road_exp(:,1)/sum(tot_dir_road_exp,2),tot_dir_road_exp(:,2)/sum(tot_dir_road_exp,2),tot_dir_road_exp(:,3)/sum(tot_dir_road_exp,2)];
tot_IEA_road=shares_dir_road*IEA_road_dir(ctry,:);


for r=1:3
shares_road_exp(:,r)=dir_road_exp(:,r)./ repmat((sum(dir_road_exp(:,r),1)),numb_hh,1);
EF_road_dir(:,r)= shares_road_exp(:,r).*repmat((tot_IEA_road(:,r)),numb_hh,1);
end


%% Concatenating direct and indirect energy footprints

Final_Direct_EF_sec_hh= [ EF_resid_dir,EF_road_dir];%5528x11
Final_Direct_EF_hh=sum(Final_Direct_EF_sec_hh,2);
%sectors that are needed to be added to the indirect energy
original_dir_sec_to_add = [55,56,58,61,62,63,131,57,59,60];% the Farm by-products sector will  be added in the end as I took only direct for it!
dir_sec_to_add = [1,2,3,4,5,6,8,9,10,11];% the 7th (Farm by-products) is taken out and will be added in the end

Final_EF_sec_per_hh=Final_indir_EF_sec_per_hh';
Final_EF_sec_per_hh(:,55)=Final_EF_sec_per_hh(:,55)+Final_Direct_EF_sec_hh(:,1);
Final_EF_sec_per_hh(:,56)=Final_EF_sec_per_hh(:,56)+Final_Direct_EF_sec_hh(:,2);
Final_EF_sec_per_hh(:,58)=Final_EF_sec_per_hh(:,58)+Final_Direct_EF_sec_hh(:,3);
Final_EF_sec_per_hh(:,61)=Final_EF_sec_per_hh(:,61)+Final_Direct_EF_sec_hh(:,4);
Final_EF_sec_per_hh(:,62)=Final_EF_sec_per_hh(:,62)+Final_Direct_EF_sec_hh(:,5);
Final_EF_sec_per_hh(:,63)=Final_EF_sec_per_hh(:,63)+Final_Direct_EF_sec_hh(:,6);
Final_EF_sec_per_hh(:,131)=Final_EF_sec_per_hh(:,131)+Final_Direct_EF_sec_hh(:,8);
Final_EF_sec_per_hh(:,57)=Final_EF_sec_per_hh(:,57)+Final_Direct_EF_sec_hh(:,9);
Final_EF_sec_per_hh(:,59)=Final_EF_sec_per_hh(:,59)+Final_Direct_EF_sec_hh(:,10);
Final_EF_sec_per_hh(:,60)=Final_EF_sec_per_hh(:,60)+Final_Direct_EF_sec_hh(:,11);

%sectors that need to replace cos they do not have indirect energy
dir_sec_to_replace=[7];

Final_EF_sec_per_hh=[Final_EF_sec_per_hh, Final_Direct_EF_sec_hh(:,7)];


%%   AGGREGATING   %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%calculating EF aggregated to 50 mix - categories
Viet_47_to_Viet_135_mask=xlsread('O:/eemba/Data/Vietnam/Viet_COICOP_to_GTAP.xlsx','VIETaggr_COICOP','C4:EG50');
Final_EF_sec_per_hh_47COICOP = Final_EF_sec_per_hh * Viet_47_to_Viet_135_mask';
 

% and now you can use the mask to get to 12 categories
Viet_12_to_Nep_135_mask=xlsread('O:/eemba/Data/Vietnam/Viet_COICOP_to_GTAP.xlsx','VIETaggr_COICOP','C54:EG65');
Final_EF_sec_per_hh_12COICOP = Final_EF_sec_per_hh * Viet_12_to_Nep_135_mask';

red_viet12_mask = Viet_12_to_Nep_135_mask';
red_viet12_mask=red_viet12_mask ((1:134),:);
FinIND_EF_sec_hh_12COICOP = Final_indir_EF_sec_per_hh'* red_viet12_mask;
%%%% Aggregating Expenditures %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%total expenditures per hh after Rasbal
exp_sec_tot = [(sum(RASbal,2)); (tot_resid_exp(:,7))]; %135x1

%expenditures before rasbal used for finding shares
Ctry_exp_incl_straw=[Ctry_exp,resid_exp(:,7)];%
for c=1:135%total number of coicop categories including straw not included in the first exp calclulation (indirect)
shares_Ctry_exp_incl_straw(:,c)=Ctry_exp_incl_straw(:,c)./repmat ((sum(Ctry_exp_incl_straw(:,c),1)),numb_hh,1);
%
Final_exp_sec_hh_tot (:,c)= repmat((exp_sec_tot(c,:)),numb_hh,1) .* shares_Ctry_exp_incl_straw(:,c);
end 

Final_EXP_sec_per_hh_47COICOP = Final_exp_sec_hh_tot  * Viet_47_to_Viet_135_mask';%
Final_EXP_sec_per_hh_12COICOP = Final_exp_sec_hh_tot  * Viet_12_to_Nep_135_mask';%5528x12


end

