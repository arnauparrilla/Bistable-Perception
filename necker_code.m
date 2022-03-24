theta = [0 1 1 0 1 0 0 0; 1 0 0 1 0 1 0 0; 1 0 0 1 0 0 1 0; 0 1 1 0 0 0 0 1; ...
         1 0 0 0 0 1 1 0; 0 1 0 0 1 0 0 1; 0 0 1 0 1 0 0 1; 0 0 0 1 0 1 1 0];

keySet = {'8','6','4','2','0'};
valueSet = {[0], [0], [0 0 0], [0 0 0], [0 0 0 0 0 0]};
P = containers.Map(keySet,valueSet);

J = 1;
b = p_brute_force(theta, P, J);
b = b.values;
disp(horzcat(b{:}));

%vector=[(:)];
%disp(vector);

%{
v1 = [1 1 1 0 0 1 0 0];
v2 = [0 1 0 1 1 1 0 0];
v3 = [1 0 0 0 1 1 0 1];

b1 = [1 1 1 0 1 0 0 0];
b2 = [1 1 0 1 0 1 0 0];
b3 = [1 1 1 0 1 0 0 0];

w1 = [1 1 1 0 0 0 0 1];
w2 = [0 1 0 1 0 1 1 0];
w3 = [0 0 1 1 0 1 1 0];

a1 = [1 0 1 0 0 1 0 1];
a2 = [1 0 0 1 1 0 0 1];
a3 = [0 1 0 1 1 0 1 0];

q1 = [1 0 0 1 0 0 0 0];
q2 = [0 1 0 0 0 0 0 1];
q3 = [0 0 1 0 0 1 0 0];

z1 = [1 0 0 0 0 0 0 1];
z2 = [0 1 0 0 0 0 1 0];
z3 = [0 0 1 0 0 1 0 0];
%}


%{
@desc function that given and integer number from 0 to 255, 
returns and 8 digit combination of 1 & -1

@param {int} num: integer from 0 to 255
@return {array} v: array with and 8 digit combination of 1 & -1
%}
function v = neckvec(num)
    vec  = dec2bin(num,8) - '0'; %8 digit binary representation
    aux_vec = 1 - vec;
    vec = vec.*-1;
    v = vec + aux_vec;
end

function k = sum_ising_prob(vec, theta) %la prob de veritat es exp això es nomes la part del sumatori (s'ha de modificar)
    k = 1/2*(vec*(theta*vec'));
end


%{
@desc compute and save into a Map obj the 14 different prob OF the necker 
cube configuration.
the cube is representet as a graphycal mod with x_i = {-1,1} and theta its
conections, we use the Issing model to get the prob ANIRIA MES A GENERAL

@param {Matrix} theta: Necker conections matrix
@param {container} P: P container with 0's arrays to store the porbs
@returns {container} Probs: P with the correct probs

%}
function Probs = p_brute_force(theta, Probs, J)
    for i = 0:255
        v = neckvec(i);
        mu = sum(v);
        k = sum_ising_prob(v,theta);
    
        switch abs(mu)
            case 8
                list = Probs('8');
                idx = 1;
                list(idx) = list(idx) + exp(J*k); %exp(J*k); ET CALCULA L'EXPONENCIAL
                Probs('8') = list;
            case 6
                list = Probs('6');
                idx = 1;
                list(idx) = list(idx) + exp(J*k);
                Probs('6') = list;
            case 4
                list = Probs('4');
                if k == 4
                    %cont = cont +1; 
                    idx = 1;
                    list(idx) = list(idx) + exp(J*k); 
                    Probs('4') = list;
                elseif k == 0
                    idx = 2;
                    list(idx) = list(idx) + 1; %la k es 0 en aquest cas 1 per visual
                    Probs('4') = list;
                end
            case 2
                list = Probs('2');
                if k == 2
                    idx = 1;
                    list(idx) = list(idx) + exp(J*k); 
                    Probs('2') = list;
                elseif k == -2
                    idx = 2;
                    list(idx) = list(idx) + exp(J*k); 
                    Probs('2') = list;
                elseif k == -6
                    idx = 3;
                    list(idx) = list(idx) + exp(J*k); 
                    Probs('2') = list;
                end
            case 0
                list = Probs('0');
                if k == 0
                    idx = 1;
                    list(idx) = list(idx) + 1; 
                    Probs('0') = list;
                elseif k == 4
                    idx = 2;
                    list(idx) = list(idx) + exp(J*k); 
                    Probs('0') = list;
                elseif k == -4
                    idx = 3;
                    list(idx) = list(idx) + exp(J*k); 
                    Probs('0') = list;
                elseif k == -12
                    idx = 4;
                    list(idx) = list(idx) + exp(J*k);%is very close to zero (matlab still save it but shows 0 
                    %when disp the Probs('0') but with "Probs('0')(4)" you get the real value)
                    Probs('0') = list;
                end
            otherwise

        end
    end
end

%{
        if mu == 8 || mu == -8
            list = Probs('8');
            idx = 1;
            list(idx) = list(idx) + J*k; %exp(J*k); ET CALCULA L'EXPONENCIAL
            Probs('8') = list;
        elseif mu == 6 || mu == -6 
            list = Probs('6');
            idx = 1;
            list(idx) = list(idx) + J*k;
            Probs('6') = list;
        elseif mu == 4 || mu == -4
            list = Probs('4');
            if k == 4
                %cont = cont +1; 
                idx = 1;
                list(idx) = list(idx) + J*k; 
                Probs('4') = list;
            elseif k == 0
                idx = 2;
                list(idx) = list(idx) + 1; %la k es 0 en aquest cas 1 per visual
                Probs('4') = list;
            end
        elseif mu == 2 || mu == -2
            list = Probs('2');
            if k == 2
                idx = 1;
                list(idx) = list(idx) + J*k; 
                Probs('2') = list;
            elseif k == -2
                idx = 2;
                list(idx) = list(idx) + J*k; 
                Probs('2') = list;
            elseif k == -6
                idx = 3;
                list(idx) = list(idx) + J*k; 
                Probs('2') = list;
            end
        elseif mu == 0 
            list = Probs('0');
            if k == 0
                idx = 1;
                list(idx) = list(idx) + 1; 
                Probs('0') = list;
            elseif k == 4
                idx = 2;
                list(idx) = list(idx) + J*k; 
                Probs('0') = list;
            elseif k == -4
                idx = 3;
                list(idx) = list(idx) + J*k; 
                Probs('0') = list;
            elseif k == -12
                idx = 4;
                list(idx) = list(idx) + J*k; 
                Probs('0') = list;
            end
        end
    end
end

%}
%Fer el GIBBS SAMPLING

%IN DEVELOPMENT
%MinorTickValues per tenir 2 tipus de x diferents   
%https://es.mathworks.com/matlabcentral/answers/22059-double-ticks-in-right-axis-plotyy
%https://es.mathworks.com/matlabcentral/answers/381959-i-would-like-to-have-two-y-axis-and-corresponding-y-axis-ticks-and-ticks-values-but-with-a-single-se

%cont
