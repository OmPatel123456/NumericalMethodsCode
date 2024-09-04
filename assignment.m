%Course : MATH∗2130
%Name : Om Patel 
%Student Number: 1228692


% begins with two initial estimates x1 and x2, then draws a line connecting
% the lines between the y values. X2 is the new estimate of the root. Use
% x1 and x2 (two most recent estimates) to approximate the roots
% iteratively 

clear
clc
close all %clear the console and close any existing output after finished 


%DEFINE VARIABLES F, X0, X1, P, NMAX

% create a symbolic variable 
syms x;

f(x)=exp(10^(-50)*x+10^(-51)) - 10^(-1)*x-1; %create the function 

rootapprox1=secant(f, 1.1*10^(-15) , 2*10^(-15),10 ,250); % one root will be at 0, so we can take test values -0.5, 0.5, to 10 digits of precision, and 
rootapprox2=secant(f, 1.1758*10^52 , 1.77*10^52 ,10,250); % call the function again in a range that included 1.1759*10^58
%note: since the numbers are so large in 

% print all the results of both calls, by printing the variables that hold
% the call results (rootapprox1, and rootapprox2)
fprintf("\n\n**final root approximation 1 : %d**\n", rootapprox1) 
fprintf("**final root approximation 2 : %d**", rootapprox2)


function [x2, f2] = secant (f, x0, x1,p, nmax)  % create the function, receiving a symbolic function f, 2 initial guesses, p for approximation, and nmax for max number of iterations 

    f0=double(f(x0)); % create a function variable f0 that holds the value of f(x0)
    f1=double(f(x1)); % create a function variable f0 that holds the value of f(x1)

    n=0; %define bound for iterations 
    re=1; %relative error set to 1, so it can decrease every time until the limit is hit (exits while loop after)

       

    while ((n<=nmax) && re>(0.5*(10^-p))) %loop through the while loop and recalculate the roots, updating x1 with the new root x2 found, and updating x0 with the previous x1. 
           
        x2= x1-f1*(x1-x0)/(f1-f0); % calculates the new root using secant methods formula 
        x0=x1; x1=x2; %update the roots 
        
        f0=f1;
        f1=double(f(x1));

        re=abs((x1-x0)/x1);
        n=n+1;  

        fprintf("\n\nroot approximation %i for x2\r: %d", n, x2) %print out the root approximations as the code progresses to keep track of the behaviour of the method and whether or not the chosen approximations are suitable 


    end 

    f2=f1; %return the updated value of f2, which holds the function 

end




Output for Secant Method: 

As seen above, the output for the secant method converges to roots that are very close to the real number. Root approximation 1 = 1e-50, which is almost 0, and root approximation 2 ends up being 1.1758*10^52 with an initial estimate of the same number. This tells us that the approximation was correct. 

ii) Newton’s Method 

%Course : MATH∗2130
%Name : Kush Shastri 
%Student Number: 1223925

% begins with one initial estimates x0

clear
clc
close all %clear the console and close any existing output after finished 


%DEFINE VARIABLES F, X0, X1, P, NMAX

% create a symbolic variable 
syms x;

f(x)=exp(10^(-50)*x+10^(-51)) - 10^(-1)*x-1; %define the function and assign it to a symbolic function f(x)

rootapprox1 =newton(f, 0, 10, 250); % one root will be at 0, so we can take test value 0 to 10 digits of precision
rootapprox2 =newton(f, 1.1759*10^52, 10, 250); % one root will be at 1.1759*10^52, so we can take test value 1.1759*10^52 to 10 digits of precision

% print all the results 
fprintf("\n\n**final root approximation 1 : %d**\n", rootapprox1) %print each final results of the two root approximations 
fprintf("**final root approximation 2 : %d**", rootapprox2)


function [xr, fr] = newton (f, x0, p, nmax)  %define newton's method. It receives a symbolic function, initial estimate of the root,
                                             % p significant digits as stopping criteria and n iterations as stopping criteria 

    xr=x0; %set the xresult to initial guess 
    n=0; %define bound for iterations 
    re= 1; %relative error 

    g=diff(f);  %differentiate the function f for newton's method 

    
    while (n<nmax && re>=double(0.5*(10^-p))) % a while loop will keep on iterating and updating the root (xr) until the relative error to the true value is less than the threshold or it has reached the max iterations           
        xrold=xr;
        xr=double(xr-f(xr)/g(xr)); %call the newton's method formula to update the root 
        re=abs((xr-xrold)/xr); %recacualte the relative error
        n=n+1;  

        fprintf("\n\nroot approximation %i for x0\r: %d", n, xr) %print each approximation as the function iterates through the roots 
    end 
    
    xr=double(xr);
    fr=double(f(xr));

end










The results for the Newton's method worked as expected, because the initial guess made was almost exactly what the function roots are supposed to be. Newton’s method only works when a root estimate is close to the real. 






iii) Fixed Point Method 

%Course : MATH∗2130
%Name : Kush Shastri 
%Student Number: 1223925


clear
clc
close all %clear the console and close any existing output after finished 


%DEFINE VARIABLES F, X0, X1, P, NMAX

% create a symbolic variable 
syms x;

f(x)=exp(10^(-50)*x+10^(-51)) - 10^(-1)*x-1; %create the function 

rootapprox1 =fixedpoint(f, 0, 10, 250); % one root will be at 0, so we can take test value 0, to 10 digits of precision 
rootapprox2 =fixedpoint(f, 1.1759*10^52, 10, 250);  % one root will be at 0, so we can take test value 1.1759*10^52 to 10 digits of precision. We will use the approximations that are as close to the actual value as possible because this makes testing easier. 

% print all the results 
fprintf("\n\n**final root approximation 1 : %d**\n", rootapprox1)
fprintf("**final root approximation 2 : %d**", rootapprox2)


%this function uses the fixed point method to take an initial value
%estimation of the root, along with a function and stopping criteria to get
%the real root. 
function [xr, fr] = fixedpoint (f, x0, p, nmax)   

%inputs: function, initial guess, p digits of approx, nmax iterations 
%outputs: 

    xr=x0;
    n=0; %define bound for iterations 
    re= 1; %relative error 

    syms x;
    g(x)=f+x;  %create the symbolic function g(x)=f+x which is what lets us do the fixed point calculations 

    
    while (n<nmax) && (re>=double(0.5*(10^-p)) ) %while our number iterations is less than the max allowed and the relative error is not to the precision that is entered (p), keep iterating
                  
        xrnext=double(g(xr)); % for fixed point method, set the next root approximation to the fixed point formula of the function + x        
        re=abs((xrnext-xr)/xrnext); %calculate the relative error 

        xr = xrnext;
        n=n+1;  %increase n by 1 so we can iterate 

        fprintf("\n\nroot approximation %i for x0\r: %d", n, xrnext)
    end 
    
    xr=double(xr);
    fr=double(f(xr));
    
end





Output: 
From this output, we can see that the fixed point method is accurate for the first root, but is inaccurate for the second root. Since the first root is at coordinate (0,0), it is a fixed point, so the method will get the root correct. According to the notes, this method calculates the roots based on finding fixed points of the function g=f(x)+x. If there are no fixed points within the interval, which is the case for root approximation 2, there will be an inaccurate result. 


Is it possible to make the “Simple” Fixed Point Method work in order to approximate all the roots of f? If not, explain when and how you were able to anticipate the failure of this strategy.

It is not possible to make a simple fixed point method work to approximate all roots of f because we can see that the second root approximation was off by a lot. The answer should have been 1.1759*10^52, yet we get approximately 5.57*10^40. The initial error needs to be very small for fixed point approximation to work. When this function g(x)=f(x)+x has a fixed point within a small range of an interval, the method works to converge towards that root. Due to the nature of this varying data, if the function does not converge properly near the roots, there may be an issue in approximating them. Fixed point method may diverge to the wrong root if there are multiple fixed points within an interval. Anticipating if this method will fail or not depends on analyzing which function we are working with. If it seems that there are multiple roots within a small range, it is important to be careful in how we approximate these roots, or we will converge to incorrect solutions. It is a good idea to consider a method with multiple initial guesses like secant method if there are a lot of roots. 






Question 2a: 

%Course : MATH∗2130
%Name : Kush Shastri 
%Student Number: 1223925


clear
clc
close all %clear the console and close any existing output after finished 


%DEFINE matrix VARIABLES Depth and Temp (p and q)
T = [0 2.3 4.9 9.1 13.9 18.3 22.9 25.1 27.2]; % x values
Z = [22.8 22.8 22.8 21.6 13.2 11.4 11.1 11.1 11.1]; %y values 


%print out the plot of the polynomial over the interval [0,27.2]
f=newton_interpolation(T,Z);
fprintf("\n\npolynomial interpolation is: ");
disp(f);

figure (1);
plotInt= [0, 27.2];
fplot(f, plotInt)
title ("Newton's Poynomial Approximation");
xlabel('Depth z (in m)') 
ylabel('Temperature t (in degrees C)')



%print out the plot of the polynomial's derivative  over the interval [0,27.2]
g=abs(diff(f));
fprintf("\n\n newton's polynomial approximation differentiated approximation is: ");
disp(g);

figure (2);
plotInt= [0, 27.2];
fplot(g, plotInt)
title ("Newton's differentiated Approximation");
xlabel('Depth z (in m)') 
ylabel('Temperature t (in degrees C)')






function [f] = newton_interpolation (T,Z)
% the function F initializes a matrix M to get data from
%this function creates a Newton polynomial and returns it 


    n=length(T); % calculates how many data points there are in the set 
    M=zeros(n);  % zero fills the matrix so that the entries can be filled using the data Z and T

    for k=1:n
        M(k,k)=Z(k); %initializes the diagonal entries, which hold the data points y1,y2,....yn
    end


    
    for i=1:n-1
        for j=1:n-i
            M(j,j+i)= (M(j+1,j+i)-M(j,j+i-1))/(T(j+i)-T(j)); % fill in the matrix M with the divided differences in the position (j,j+i)
        end

    end


    %------------------------------------------------------------------------------------------
    %code the newton's polynomial function 

    syms x; %symbolic variable x used for the polynomial 
    syms f(x); %symbolic function that will hold the final polynomial 
    syms newprod(x); %newprodx is the product, updated after each iteration
    newprod(x)=1; %initialiiy set the product to 1 so 
    f(x)=M(1,1); %initializing value for the final polynomial 
  


    for i= 2:n
        %constructs the newton polynoimal with coefficients term-by-term 
        newprod(x)= newprod(x)*(x-T(i-1));
        f(x)= f(x) + M(1,i)*newprod(x);
        % the coefficients start at 1,  as we want the ones off the top of
        % the tree 

    end

   %------------------------------------------------------------------------------------------

end 

