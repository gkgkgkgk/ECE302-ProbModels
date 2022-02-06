
%% Q1
clear;
clc;
numSimulations = 1000000;

% a)
% P(X = 18) = 1/6 * 1/6 * 1/6 = 4.63e-3
disp("1A:");
calculatedProb = (1/6) * (1/6) * (1/6);
rolls = 1:numSimulations;

for i = 1:numSimulations
    rolls(i) = diceRoll(6) + diceRoll(6) + diceRoll(6);
end

p = sum(rolls == 18) / numSimulations;
disp("Calculated Probability of an 18: " + calculatedProb);
disp("Simulated Probability of an 18: " + p);


% b)
% P(Xfun = 18) = 1-(215/216)^3 = 1.38e-2
disp("1B:");

calculatedProb = 1-(215/216)^3;

for i = 1:numSimulations
    r(1) = diceRoll(6) + diceRoll(6) + diceRoll(6);
    r(2) = diceRoll(6) + diceRoll(6) + diceRoll(6);
    r(3) = diceRoll(6) + diceRoll(6) + diceRoll(6);

    rolls(i) = max(r);
end

p = sum(rolls == 18) / numSimulations;
disp("Calculated Probability of an 18 using The Fun Method: " + calculatedProb);
disp("Simulated Probability of an 18 using The Fun Method: " + p);

% c)
% P(Xn = 18), 1 <= n <= 6 = P(Xfun)^6 = 6.985e-12
disp("1C:");

calculatedProb = (1-(215/216)^3)^6;

c = zeros(1, numSimulations);

for i = 1:numSimulations
    for j = 1:6
        r(1) = diceRoll(6) + diceRoll(6) + diceRoll(6);
        r(2) = diceRoll(6) + diceRoll(6) + diceRoll(6);
        r(3) = diceRoll(6) + diceRoll(6) + diceRoll(6);
        c(i) = c(i) + max(r);
    end
end

p = sum(c == 108) / numSimulations;

disp("Calculated Probability of an 18 using The Fun Method: " + calculatedProb);
disp("Simulated Probability of an 18 using The Fun Method: " + p);

% d)
% P(Xn = 9), 1 <= n <= 6 = P(Xfun = 9)^6 = (25/216)^6
disp("1D:");

calculatedProb = (25/216)^6;
c = zeros(6, numSimulations);

for i = 1:numSimulations
    for j = 1:6
        r(1) = diceRoll(6) + diceRoll(6) + diceRoll(6);
        r(2) = diceRoll(6) + diceRoll(6) + diceRoll(6);
        r(3) = diceRoll(6) + diceRoll(6) + diceRoll(6);
        c(j, i) = max(r);
    end
end

p = 0;

for i = 1:numSimulations
    if c(:, i) == [9 9 9 9 9 9]
        p = p + 1;
    end
end

p = p / numSimulations;

disp("Calculated Probability of an 18 using The Fun Method: " + calculatedProb);
disp("Simulated Probability of an 18 using The Fun Method: " + p);

%% 2
clear;
clc;
close all;

numSimulations = 10000;

% a)
disp("2A:");

trolls = 1:numSimulations;

for i = 1:numSimulations
    trolls(i) = diceRoll(4);
end

a = mean(trolls); % haha mean trolls. 
avgTrollHP = (1/4) + (2/4) + (3/4) + (4/4);

disp("Calculated AVG hitpoints per troll: " + avgTrollHP);
disp("Simulated AVG hitpoints per troll: " + a);

% Likelihood to do more than 3 points of damage:
% P(X > 3) = P(X = 4) = 1/4

p = (1/4);
spells = zeros(1, numSimulations);

for i = 1:numSimulations
    spells(i) = diceRoll(2) + diceRoll(2);
end

c = sum(spells == 4) / numSimulations;

disp("Calculated likelihood to do more than 3 points of damage: " + p);
disp("Simulated likelihood to do more than 3 points of damage: " + c);

% b)
trollsPMF(1) = sum(trolls == 1)/numSimulations;
trollsPMF(2) = sum(trolls == 2)/numSimulations;
trollsPMF(3) = sum(trolls == 3)/numSimulations;
trollsPMF(4) = sum(trolls == 4)/numSimulations;

spellsPMF(1) = sum(spells == 2)/numSimulations;
spellsPMF(2) = sum(spells == 3)/numSimulations;
spellsPMF(3) = sum(spells == 4)/numSimulations;

figure('Name', "2B");
subplot(2, 1, 1);
stem(trollsPMF);
xticks([1, 2, 3, 4]);
title("PMF for HP of Trolls")
xlabel("Troll HP");
ylabel("Probability");
subplot(2, 1, 2);

stem(spellsPMF);
xticks([2, 3, 4]);
title("PMF for FIREBALL Damage Output")
xlabel("FIREBALL Damage");
ylabel("Probability");

% c)
% P(Y > Xn), 1 <= n <= 6 where Y is the damage value and X is the troll hp.
% P(Y > Xn), 1 <= n <= 6 = (1/4)(2/4)^6 + (1/2)(3/4)^6 + (1/4)(1)^6

disp("2C:");

survived = zeros(1, numSimulations);
p = (1/4)*(2/4)^6 + (1/2)*(3/4)^6 + (1/4)*(1)^6;

for i = 1:numSimulations
    fireball = diceRoll(2) + diceRoll(2);
 
    for j = 1:6
        hp = diceRoll(4);    
        if hp > fireball
            survived(i) = survived(i) + 1;
        end
    end
end

survivorCount = sum(survived == 0) / numSimulations;

disp("Calculated Probability that no trolls survived:" + p);
disp("Simulated Probability that no trolls survived:" + survivorCount);

% d)
disp("2D:");

% Scenarios for only one troll surviving:
%
% Dice roll of 3:
% Troll Survivor HP = 4, all other trolls are 3 or below
%
% Dice roll of 2:
% Troll Survivor HP = 4, all other trolls are 2 or below
%
% Dice roll of 2:
% Troll Survivor HP = 3, all other trolls are 2 or below

% total probablity:
p = (1/4)*(1/4)*(3/4)^5 + (1/4)*(1/4)*(2/4)^5 + (1/4)*(1/4)*(2/4)^5;
% probability when troll hp is 4:
hp4 = ((1/4)*(1/4)*(3/4)^5 + (1/4)*(1/4)*(2/4)^5) / p;
% probability when troll hp is 3:
hp3 = ((1/4)*(1/4)*(2/4)^5) / p;

expectedHP = (4 * hp4 + 3 * hp3);

survivor = [];

for i = 1:numSimulations
    fireball = diceRoll(2) + diceRoll(2);
    
    troll = 0;
    
    for j = 1:6
        hp = diceRoll(4);
        if hp > fireball
            if troll == 0
                troll = hp;
            else
                troll = 0;
                break;
            end
        end
    end
    
    if troll > 0
        survivor(end + 1) = troll;
    end

end

disp("Calculated HP of remaining troll:" + expectedHP);
disp("Simulated HP of remaining troll:" + mean(survivor));

% e)
% E[X] = (10/20)(E[2d6]) + (10/40)(E[1d4]) = (1/2)(7) + (1/4)(2.5) =  4.125
damage = zeros(1, numSimulations);
expectedDMG = (1/2)*(7) + (1/4)*(2.5);

for i = 1:numSimulations
    if(randi(20) >= 11) % attack with sword
        damage(i) = damage(i) + diceRoll(6) + diceRoll(6);
        
        if(randi(20) >= 11)
            damage(i) = damage(i) + diceRoll(4);
        end
    end
end

disp("Calculated damage from Shedjam:" + expectedDMG);
disp("Simulated damage from Shedjam:" + mean(damage));

% f) 
% 1. awesome assignment

function d = diceRoll(sides)
    d = randi(sides);
end