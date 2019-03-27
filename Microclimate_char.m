%Importo i dati grezzi di T, RH e tempo
T = [20.0;23.1;23.1;23.2;23.1;23.1;23.1;23.2;23.1;23.1;23.1;23.2;23.1;23.1;23.1;23.2];
RH = [45;51;51;51;51;55;56;57;60;60;60;60;60;60;60;65];
t = [0:15:225]';
n = size(T,1);

%Converto in EMC: EMC è una matrice che ha il tempo come prima colonna e l'EMC come seconda

EMC = zeros(n,2);
for j = 1:n
    EMC(j,2) = 100*(-0.0991*tanh(-10.0529*RH(j)/100+0.0034*T(j)+9.9448)+16.0049*tanh(0.5071*RH(j)/100-0.0014*T(j)-2.8132)+16.0152);
    EMC(j,1) = t(j);
end

%Eseguo la media mobile ogni npti_media punti

npti_media = input('Ogni quanti punti vuoi mediare?');

if mod(n,npti_media) == 0
    continue
else
    EMC = EMC(1:end-mod(n,npti_media),:);
end

EMC_smooth = zeros(n/npti_media,2);
m = 1;
l = floor(npti_media/2)+1;

for j = 1:(n/npti_media)
    %Tempi
    EMC_smooth(j,1) = EMC(l,1);
    %EMC
    EMC_smooth(j,2) = mean(EMC(m:(m+npti_media-1),2));
    
    m = m + npti_media;
    l = l + npti_media;  
end

%Inserisco i due parametri per l'analisi del microclima
delta = input('Inserire il valore di delta. \n\delta = ');
Tmin = input('Inserire il valore di Tmin (in ore). \n\Tmin = ');

%Costruisco i PERIODI successivi, sulla base della sola scelta di delta

periods = zeros(1,5); %Prima colonna: EMC media, Seconda colonna: t iniz., Terza colonna: Delta t, Quarta colonna: EMC_min, Quinta colonna: EMC_max

iniz = 1;
numper = 1;
while iniz<size(EMC_smooth,1)
    for j = iniz:size(EMC_smooth,1)
        ev1 = abs(EMC_smooth(j+1,2)-EMC_smooth(j,2));
        ev2 = abs(max(EMC_smooth(iniz:j+1,2))-mean(EMC_smooth(iniz:j+1,2))); 
        ev3 = abs(min(EMC_smooth(iniz:j+1,2))-mean(EMC_smooth(iniz:j+1,2)));
        
        if (ev1<delta && ev2<delta && ev3<delta)
            continue
        else
            periods(numper,1) = mean(EMC_smooth(iniz:j,2));
            periods(numper,2) = EMC_smooth(iniz,1);
            periods(numper,3) = EMC_smooth(j,1)-EMC_smooth(iniz,1);
            periods(numper,4) = min(EMC_smooth(iniz:j,2));
            periods(numper,5) = max(EMC_smooth(iniz:j,2));
            periods = [periods;0 0 0 0 0];
            iniz = j+1;
            numper = numper+1;
            break
        end
        if j == size(EMC_smooth,1)
            periods(numper,1) = mean(EMC_smooth(iniz:j,2));
            periods(numper,2) = EMC_smooth(iniz,1);
            periods(numper,3) = EMC_smooth(j,1)-EMC_smooth(iniz,1);
            periods(numper,4) = min(EMC_smooth(iniz:j,2));
            periods(numper,5) = max(EMC_smooth(iniz:j,2));
            iniz = j;
        else
            continue
        end
    end
end
if periods(end,:) == 0
    periods = periods(1:end-1,:);
end

%Tra i PERIODI costruiti, distinguo quali sono PLATEAUX e quali appartengono invece a PERIODI DI TRANSIZIONE. La matrice disc_mat identifica con valore True i plateaux, con valore False le transizioni


disc_mat = false(size(periods,1),1);

for j = 1:size(periods,1)
    if periods(j,3) >= Tmin
        disc_mat(j) = 1;
    end
end

plateaux = periods(disc_mat,:);
trans_mat = periods(!disc_mat,:);

%Unisco le transizioni consecutive a formare i periodi di transizione finali

if size(trans_mat,1) > 1
    transitions = zeros(1,4); %Prima colonna: DELTA_EMC, Seconda colonna: t iniz., Terza colonna: Delta_t, Quarta colonna: p

    numtransition = 1;
    iniz = 1;
    while iniz < size(trans_mat,1)
        for j = iniz:size(trans_mat,1)
            if j == size(trans_mat,1)
                transitions(numtransition,1) = trans_mat(j,5)-trans_mat(iniz,4);
                transitions(numtransition,2) = trans_mat(iniz,2);
                transitions(numtransition,3) = trans_mat(j,2)+trans_mat(j,3)-trans_mat(iniz,2);
                transitions(numtransition,4) = transitions(numtransition,1)/transitions(numtransition,3);
                iniz = j;
                break
            else
                continue
            end 
            if (trans_mat(j+1,2)-trans_mat(j,2)) <= Tmin
                continue
            else
                transitions(numtransition,1) = trans_mat(j,5)-trans_mat(iniz,4);
                transitions(numtransition,2) = trans_mat(iniz,2);
                transitions(numtransition,3) = trans_mat(j,2)+trans_mat(j,3)-trans_mat(iniz,2);
                transitions(numtransition,4) = transitions(numtransition,1)/transitions(numtransition,3);
                transitions = [transitions;0 0 0 0];
                iniz = j+1;
                numtransition = numtransition + 1;
                break
            end
            if j == size(trans_mat,1)
                transitions(numtransition,1) = trans_mat(j,5)-trans_mat(iniz,4);
                transitions(numtransition,2) = trans_mat(iniz,2);
                transitions(numtransition,3) = trans_mat(j,2)+trans_mat(j,3)-trans_mat(iniz,2);
                transitions(numtransition,4) = transitions(numtransition,1)/transitions(numtransition,3);
                iniz = j;
            else
                continue
            end
        end
    end
    if transitions(end,:) == 0
        transitions = transitions(1:end-1,:);
    end
elseif size(trans_mat,1) == 1
    transitions = zeros(1,4);
    transitions(1,1) = trans_mat(1,5)-trans_mat(1,4);
    transitions(1,2) = trans_mat(1,2);
    transitions(1,3) = trans_mat(1,3);
    transitions(1,4) = transitions(1,1)/transitions(1,3);
else
    transitions = [];
end