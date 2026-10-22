function [cea_output]=wrp_rkt_AP_HTPB_tran(p,pi_pe,temp_grain, APfrac, HTPBfrac)

%-----------------------INPUT------------------------
% p -> pressione in camera 
% aeat -> rapporto exp ugello 
% temp_grain -> temperatura grano
%----------------------OUTPUT------------------------
%l'output è una struct: 

% assegno a variabile l'output di un wrapper x=wrp_*tipoDiwrap*(p,ae/at, temp_grain) 
% x.output -> raggiongo i campi output 
% x.output.eql -> raggiongo i campi per soluzione full shifting
% x.output.frz -> raggiongo i campi per soluzione shifting+frozen 
% x.output.eql.field -> l'output segue la natura degli output CEA , su un
% problema rocket x.output.eql.fieldName restituisce [field_in_chamber , field_in_throat, field_in_exit]

%VARIABILI CHE NON DIPENDONO DA SHIFTING O FROZEN (ex o/f , phi, %fuel
%, Requiv) SONO RAGGIUNGIBILI DA x.output.fieldName 
% ex o/f=x.output.oxfl 

%Gli output disponibili sono molteplici , consiglio di eseguire
%x=wrp_simpleCase() e esplorare la struttura della soluzione dal riquadro
%workspace
%-------------------------------------------------------------------------



oxFracFuel= APfrac/HTPBfrac; 
wtHTPB= 100; 


cea_output = CEA('problem', 'rocket', 'frozen', 'o/f', oxFracFuel, ...
    'case', 'BASELINE_RKT', 'p,bar', p, 'Pi/Pe', pi_pe , ...
    'reactants', ...
    'fuel','HTPB ', 't,k',temp_grain,'h,kj/mol',-58 , 'C' ,7.075 ,'H' ,10.65 ,'O' ,0.223 ,'N' ,0.063,'wt%', wtHTPB, ...
    'oxid', 'NH4CLO4(I)', 't,k', temp_grain, 'wt%', 100, ...
    'output' , 'massf' ,'tran','end'  ); 


end 
