%% Compare Execution Times
TimeGen=load("TimeGen");
TimeRan=load("TimeRan");
TimeBB=load("TimeBB");
TimeBBMin=load("TimeBBMin");

figure(1);
subplot(3, 1,1);
plot(TimeGen(:,1), log10(TimeGen(:,2)), 'b','LineWidth',2);
title('Execution Time of Brute Force General Permutations');
xlabel("N");
ylabel("Execution Time log10(sec)");

subplot(3, 1,2);
plot(TimeRan(:,1), log10(TimeRan(:,2)), 'b','LineWidth',2);
title("Execution Time of Brute Force Random Permutations");
xlabel("N");
ylabel("Execution Time log10(sec)");

subplot(3, 1,3);
plot(TimeBB(:,1), log10(TimeBB(:,2)), 'b','LineWidth',2);
title("Execution Time of Branch-and-Bound");
xlabel("N");
ylabel("Execution Time log10(sec)");


figure(2);
hold on;
title("Execution Time of Branch-and-Bound");
plot(TimeBBMin(:,1), log10(TimeBBMin(:,2)), 'g','LineWidth',2);
plot(TimeBB(:,1), log10(TimeBB(:,2)), 'b','LineWidth',2);
legend('Branch-and-Bound Before', 'Branch-and-Bound Improved');
xlabel("N");
ylabel("Execution Time log10(sec)");
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% For N=13 and N=14 graphics

GenN13=load("GenN13");
RanN13=load("RanN13");

GenN14=load("GenN14");
RanN14=load("RanN14");

figure(3);
c=GenN13(:,1);
n_c=GenN13(:,2);
med=sum(n_c.*c)/sum(n_c);
var=sum(n_c.*(c-med).^2)/sum(n_c);
gaus=(1/sqrt(2*pi*var))*exp(-(c-med).^2/(2*var))*sum(n_c);
hold on;
bar(GenN13(:,1), GenN13(:,2), 'b');
plot(c, gaus, 'k','LineWidth',1.5);
legend("Gaussian", "Histogram");
hold off;title('Costs of Brute Force General Permutations for N=13');
xlabel("Costs");
ylabel("Number of Times x Appears");

figure(4)
c=RanN13(:,1);
n_c=RanN13(:,2);
med=sum(n_c.*c)/sum(n_c);
var=sum(n_c.*(c-med).^2)/sum(n_c);
gaus=(1/sqrt(2*pi*var))*exp(-(c-med).^2/(2*var))*sum(n_c);
hold on;
bar(RanN13(:,1), RanN13(:,2), 'g');
plot(c, gaus, 'k','LineWidth',1.5);
legend("Gaussian", "Histogram");
hold off;
title('Costs of Brute Force Using Random Permutations for N=13');
xlabel("Costs");
ylabel("Number of Times x Appears");

figure(5);
c=GenN14(:,1);
n_c=GenN14(:,2);
med=sum(n_c.*c)/sum(n_c);
var=sum(n_c.*(c-med).^2)/sum(n_c);
gaus=(1/sqrt(2*pi*var))*exp(-(c-med).^2/(2*var))*sum(n_c);
hold on;
bar(GenN14(:,1), GenN14(:,2), 'b');
plot(c, gaus, 'k','LineWidth',1.5);
legend("Gaussian", "Histogram");
hold off;
title('Costs of Brute Force General Permutations for N=14');
xlabel("Costs");
ylabel("Number of Times x Appears");

figure(6);
c=RanN14(:,1);
n_c=RanN14(:,2);
med=sum(n_c.*c)/sum(n_c);
var=sum(n_c.*(c-med).^2)/sum(n_c);
gaus=(1/sqrt(2*pi*var))*exp(-(c-med).^2/(2*var))*sum(n_c);
hold on;
bar(RanN14(:,1), RanN14(:,2), 'g');
plot(c, gaus, 'k','LineWidth',1.5);
legend("Gaussian", "Histogram");
hold off;
title('Costs of Brute Force Using Random Permutations for N=14');
xlabel("Costs");
ylabel("Number of Times x Appears");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Minimums and maximums

MaxCostBB=load("MaxCostBB");
MaxCostRan=load("MaxCostRan");
MaxCostGen=load("MaxCostGen");
MinCostGen=load("MinCostGen");
MinCostRan=load("MinCostRan");
MinCostBB=load("MinCostBB");

figure(7);
hold on;
plot(MinCostGen(:,1), MinCostGen(:,2), 'b','LineWidth',3);
plot(MinCostRan(:,1), MinCostRan(:,2), 'g','LineWidth',1.5);
plot(MinCostBB(:,1), MinCostBB(:,2), 'k','LineWidth',1.5);
legend("Geneneral Permutations", "Random Permutations", "Branch-and-Bound");
hold off;
title("Minimum Costs for N's");
ylabel("Costs");
xlabel("N");

figure(8);
hold on;
plot(MaxCostGen(:,1), MaxCostGen(:,2), 'b','LineWidth',3);
plot(MaxCostRan(:,1), MaxCostRan(:,2), 'g','LineWidth',1.5);
plot(MaxCostBB(:,1), MaxCostBB(:,2), 'k','LineWidth',1.5);
legend("Geneneral Permutations", "Random Permutations", "Branch-and-Bound");
hold off;
title("Maximum Costs for N's");
ylabel("Costs");
xlabel("N");