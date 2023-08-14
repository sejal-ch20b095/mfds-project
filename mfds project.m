%QUESTION 1
%LINEAR ALGEBRA
%% 
clear all
for iter=1:5
    fprintf('\n');
    fprintf('\nQ1V%d\n', iter);
    fprintf("\n")
    x1=randi([1 6], 3, 1);
    x2=randi([1 6], 3, 1);
    var1=randi([1, 6]);
    var2=randi([-5, 5]);
    while var2==0
        var2=randi([-5, 5]);
    end
    x3=var1*x1 + var2*x2;
    X=[x1 x2 x3];
    rank(X);
    %add error  
    %
    std(1)=0.15;
    std(2)=0.15;
    std(3)=0.15;
    %
    L=diag(std);
    for j = 1:3
        error = randn(3,1);
        Xmeas(j,:) = X(j,:) + (L*error)';
    end
    rank(Xmeas);
    [v, lamb]=eig(X'*X);
    relation=-1*v(:, 1)/v(3, 1);
    weights=randi([2, 5], 2, 1);
    %display the question
    fprintf("matrix A(3*3) represents data from a process, where each column is a variable. It is known that there is a single relation between the variables \n")
    fprintf("[")
    fprintf("\n[%.2f   %.2f   %.2f]", Xmeas)
    fprintf("]")
    fprintf("\nColumns 1, 2, 3 are represented as v1, v2, v3 respectively. The linear relation between them is given by ")
    fprintf("\na*v1 + b*v2 -v3=0")
    fprintf("\nfind %.0fa + %.0fb\n", weights(1), weights(2))
    %generate the options
    ans=weights(1)*relation(1)+weights(2)*relation(2);
    format bank
    ops=[ans ans-1 ans+1 ans-2];
    choices=ops;
    choices_rand=choices(randperm(4));
    for i=1:4
        fprintf("\n %c.", i-1+'a');
        fprintf("%.2f", choices_rand(i));
        %disp(choices_rand(i));
        if(choices_rand(i)==choices(1))
            Ansopt=char(i-1+'a');
        end
    end
    fprintf("\n")
    fprintf('\n correct option is: option %c \n', Ansopt);
    %SOLUTION
    fprintf("SOLUTION\n")
    fprintf("Concept used:\n")
    fprintf("The eigenvalues of the matrix (X^T *X) gives the amount of variance captured by each variable.\n")
    fprintf("The eigen vector corresponding to the smallest eigen value gives the required single linear relation\n")
    fprintf("The eigen vector corresponding to the lowest eigen value in this question is [%.2f %.2f %.2f]' \n", v(1, 1), v(2,1), v(1, 3))
    fprintf("normalizing with the coefficient corresponding to the third variable we get the value of a and b to be %.0f %.0f \n", relation(1), relation(2))
end


%% QUESTION 2
%LINEAR ALGEBRA
%%
clear all;
for iter=1:5
    fprintf('\n')
    fprintf('\nQ2V%d\n ',iter);
    ncol=3;
    A=zeros(3, 3);
    A(1, 1)=randi([-5, 5]);
    A(1, 2)=randi([-5, 5]);A(2, 1)=A(1, 2);
    A(2, 3)=randi([-5, 5]);A(3, 2)=A(2, 3);
    A(3, 1)=randi([-5, 5]);A(1, 3)=A(3, 1);
    A(2, 2)=randi([-5, 5]);
    A(3, 3)=randi([-5, 5]);
    %A=randi([-5, 5], ncol, ncol);
    determ=det(A);
    [eig_vec, eig_val] = eig(A);
    det_eig = det(eig_vec);
    while determ==0 | det_eig == 0
        A=zeros(3, 3);
        A(1, 1)=randi([-5, 5]);
        A(1, 2)=randi([-5, 5]);A(2, 1)=A(1, 2);
        A(2, 3)=randi([-5, 5]);A(3, 2)=A(2, 3);
        A(3, 1)=randi([-5, 5]);A(1, 3)=A(3, 1);
        A(2, 2)=randi([-5, 5]);
        A(3, 3)=randi([-5, 5]);
        determ=det(A);
        [eig_vec, eig_val] = eig(A);
        det_eig = det(eig_vec);
        if(determ~=0 && det_eig~=0)
            break;
        end
    end
    A = eig_vec*eig_val*inv(eig_vec);
    a = randi([-1, 1]);
    b = randi([-1, 1]);
    c = randi([-1, 1]);
    d = randi([-1, 1]);
    %creating matrix C
    C = [a*A b*A];
    C = [C;[c*A,d*A]];
    r = rank(C);
    B = r*A;
    %Printing the question
    fprintf("Consider a matrix A with the following eigen values \n")
    fprintf("[%.2f %.2f %.2f]'\n", eig_val(1, 1), eig_val(2, 1), eig_val(3, 1))
    fprintf("and corresponding eigenvectors \n")
    fprintf("[%.2f %.2f %.2f]'\n", eig_vec(1, 1), eig_vec(2, 1), eig_vec(3, 1))
    fprintf("[%.2f %.2f %.2f]'\n", eig_vec(1, 2), eig_vec(2, 2), eig_vec(3, 2))
    fprintf("[%.2f %.2f %.2f]'\n", eig_vec(1, 3), eig_vec(2, 3), eig_vec(3, 3))
    %eig_vec(:,1)
    %eig_vec(:,2)
    %eig_vec(:,3)
    disp('respectively. Let C be a matrix defined as')
    fprintf("[[aA bA]\n")
    fprintf(" [cA dA]]")
    fprintf("where a=%.0f b=%.0f c=%.0f d=%.0f\n", a, b, c, d)
    disp('the rank of matrix C is denoted by r and matrix B is defined as B = r*A report the trace of matrix B')
    ans=B(1, 1)+B(2,2)+B(3, 3);
    fprintf('\n')
    ops=[ans ans+1 ans-1 ans+2];
    choices=ops;
    choices_rand=choices(randperm(4));
    for i=1:4
        fprintf("\n %c.", i-1+'a');
        fprintf("%.2f", choices_rand(i));
        %disp(choices_rand(i));
        if(choices_rand(i)==choices(1))
            Ansopt=char(i-1+'a');
        end
    end
    fprintf("\n")
    fprintf('\n correct option is: option %c \n', Ansopt);
    %SOLUTION
    fprintf("\nSOLUTION")
    fprintf("\nConcepts used: eigen value decomposition and rank")
    fprintf("\nmatrix A=Q*lambda*Q^-1")
    fprintf("\nFrom the given values of Q and lambda(eigen values), A can be computed")
    fprintf("\nThe matrix C has rank same as matrix A if (a=b, c=d) or (a=c, b=d)")
    fprintf("\nIn all other cases, the matrix C has twice the rank of matrix A. In this way the rank can be computed")
    fprintf("\nMatrix C's trace can be calculated as (trace(A))*rank(C)\n")
end


%% QUESION 3 OPTIMISATION
clear all;
%%
for i=1:5
    fprintf('\n');
    fprintf('\nQ3V%d\n ',i);
    fprintf('\n');
    a = randi([-10,10]);
    b = randi([-10,10]);
    c = randi([-10,10]);
    d = randi([-10,10]);
    e = randi([-10,10]);
    while a == b == c == d == e == 0
        a = randi([-10,10]);
        b = randi([-10,10]);
        c = randi([-10,10]);
        d = randi([-10,10]);
        e = randi([-10,10]);
    end

    fun = @(x) a*x(1)^3 + b*x(2)^3 - c*x(1)*x(2) + d*x(1) + e*x(2);
    nconstr = randi([1,3]);
    A = randi([-10,10],nconstr,2);
    b = randi([-100,100],nconstr,1);
    neqconstr = 1;
    if neqconstr == 1
        Aeq = randi([-10,10],1,2);
        beq = randi([5,10]);
    else
        Aeq = [];
        beq = [];
    end
    lb = rand(2,1);
    ub = rand([-10,10]);
    x0 = rand(2,1);
    while Aeq(1)/lb(1) == Aeq(2)/lb(2)
        lb = rand(2,1);
        ub = randi([5,20]);
        x0 = rand(2,1);
    end
    options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'off');
    warning('off', 'all')
    [x, fval, exitflag, output, lambda] = fmincon(fun, x0, A, b, Aeq, beq, lb, ub, @nonlcon, options);
    while fval <=0
        a = randi([-10,10]);
        b = randi([-10,10]);
        c = randi([-10,10]);
        d = randi([-10,10]);
        e = randi([-10,10]);
        while (a == b) & (b == c) &(c==d) & (d == e) & (e == 0)
            a = randi([-10,10]);
            b = randi([-10,10]);
            c = randi([-10,10]);
            d = randi([-10,10]);
            e = randi([-10,10]);
        end
        fun = @(x) a*x(1)^3 + b*x(2)^3 - c*x(1)*x(2) + d*x(1) + e*x(2);
        nconstr = randi([1,3]);
        A = randi([-10,10],nconstr,2);
        b = randi([-100,100],nconstr,1);
        neqconstr = 1;
        Aeq = randi([-10,10],1,2);
        beq = randi([5,10]);
        lb = rand(2,1);
        ub = rand([-10,10]);
        x0 = rand(2,1);
        while Aeq(1)/lb(1) == Aeq(2)/lb(2)
            lb = rand(2,1);
            ub = randi([5,20]);
            x0 = rand(2,1);
        end
        options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'off');
        warning('off', 'all')
        [x, fval, exitflag, output, lambda] = fmincon(fun, x0, A, b, Aeq, beq, lb, ub, @nonlcon, options);
    end   
    ansi = fval;
    ops = [fval randi([1,3])*fval + randi([10,20]) randi([2,5])*fval + randi([14,24]) rand*fval + randi([0,5])];
    ops = ops(randperm(length(ops)));
    for i=1:4
        fprintf("\n %c.", i-1+'a');
        fprintf("%.2f", ops(i));
        %disp(choices_rand(i));
        if ops(i)== ansi
            Ansopt=char(i-1+'a');
        end
    end
    fprintf('A company wants to minimize the cost of producing two products, P1 and P2, subject to three constraints: a limited amount of raw material, a minimum production volume, and a maximum production volume.\n')
    fprintf('The cost of producing each product is a function of the production volumes and is given by C(P1, P2) = %.2f P1^3 + %.2f P2^3 - %.2f P1P2 + %.2f P1 + %.2f P2. \n', a,b,c,d,e )
    fprintf('The raw material constraint is %.2f P1 + %.2f P2 <= %.2f, \n' , Aeq(1), Aeq(2), beq )
    fprintf('and the maximum production volume constraint is %.2f P1 + %.2f P2 <= %.2f. \n', lb(1), lb(2), ub )
    fprintf('You can choose the initial point to be (%.2f, %.2f) \n' , x0(1), x0(2))
    fprintf("\n")
    fprintf('\n correct option is: option %c \n', Ansopt);
    fprintf('Optimal solution: P1 = %.2f, P2 = %.2f\n', x(1), x(2));
    fprintf('Minimum cost: $%.2f\n', fval);

    %%Solution 
    fprintf("To solve the optimization problem using the KKT conditions, we first write the Lagrangian function: \n")
    fprintf("L(P1, P2, λ1, λ2, λ3) = P1^3 + 3P2^3 - 4P1P2 + 2P1 + 5P2 - λ1(2P1 + P2 - 100) - λ2(P1 + P2 - 20) - λ3(P1 + P2 - 50) \n")
    fprintf("The KKT conditions are: \n")
    fprintf("1. Stationarity condition: ∇L(P1, P2, λ1, λ2, λ3) = 0 \n")
    fprintf("2. Primal feasibility conditions: g1(P1, P2) = 2P1 + P2 - 100 ≤ 0, g2(P1, P2) = P1 + P2 - 20 ≥ 0, g3(P1, P2) = P1 + P2 - 50 ≤ 0 \n")
    fprintf("3. Dual feasibility conditions: λ1 ≥ 0, λ2 ≥ 0, λ3 ≥ 0 \n")
    fprintf("4. Complementary slackness conditions: λ1g1(P1, P2) = 0, λ2g2(P1, P2) = 0, λ3g3(P1, P2) = 0 \n")
    fprintf("We use the fmincon function to perform this KKT operation in MATLAB \n")
    fprintf("Optimal solution: P1 = 3.9088, P2 = 16.0912 \n")
    fprintf("Minimum cost: $12395.7592\n")



end



%% QUESTION 4 OPTIMISATION
clear all

%%
for i=1:5
    %Title
    fprintf('\n');
    fprintf('\nQ4V%d\n', i);

    %objective function
    %max height function
    par1=0.1*randi([1, 5]);
    par2=round(randi([50, 150]), -1);
    par3=round(randi([50, 150]), -1);
    par4=0.01*randi([10, 50]);

    optimfun=@(x) -1*(-1*par1*x(1)^2 +par2*x(1) +par3*x(2)-par4*x(2)^2);
    %total cost constraint
    cperweight=randi([1, 5]);
    cpervol=randi([1, 5]);
    cost=round(randi([1000, 4000]), -3);
    x0=[1, 1];
    lb=[1, 1];
    ub=[];
    A=[cpervol, cperweight];
    b=[cost];
    Aeq=[];
    beq=[];
    x0=10*lb;
    options = optimset('Display', 'off');
    [x, fval]= fmincon(optimfun,x0,A,b,Aeq,beq, lb, ub, @nonlcon1, options);
    diary off;
    fprintf("optimum combination of weight and volume is %.1f %.1f \n", x(1), x(2))
    fprintf("function value at optimum is %.1f \n", -1*fval)
    %display the question
    disp("A balloon for measuring meteorological data is sent from ground. The Scientists need to optimise the maximum height reached by altering the volume and weight of the balloon\n")
    fprintf("The maximum height as a function of height and volume is given by: h(w, v)= -%.1fw^2 + %.0fw +%.0fv -%.1fv^2 \n", par1, par2, par3, par4)
    fprintf("The minimum weight and volume must be 1. The total cost of the project is %.0f and is a function of weight and volume. \n", cost)
    fprintf("Cost per unit weight is %.0f and the cost per unit volume is %0.f \n", cperweight, cpervol)
    fprintf("Find the maximum height the balloon can go with the current constraints")
    %display options
    format short
    max_h=-1*fval;
    max_h=round(max_h, 1);
    ops=[max_h; 2*max_h; 0.75*max_h; 0.5*max_h];
    choices=ops;
    choices_rand=choices(randperm(4));
    for i=1:4
        fprintf("\n %c.", i-1+'a');
        fprintf("%.1f", choices_rand(i));
        %disp(choices_rand(i));
        if(choices_rand(i)==choices(1))
            Ansopt=char(i-1+'a');
        end
    end
    fprintf('\n correct option is: option %c \n', Ansopt); 
    disp('SOLUTION Question 3:')
    fprintf(" The optimisation function is h(w, v)= -%.1fw^2 + %.0fw +%.0fv -%.1fv^2 \n", par1, par2, par3, par4)
    fprintf("The constraints for minimum weight and volume are 1-w<=0 and 1-v<=0 \n")
    fprintf("The cost constraint is w*%.0f + v*%.0f - %.0f <= 0 \n", cperweight, cpervol, cost)
    fprintf("Thus we have 3 inequality constraint which can be solved using kkt \n")
    fprintf("solving we get the value of weight and volume to be:w=%.1f v=%.1f \n", x(1), x(2))
    fprintf("The corresponding function value is the required answer h=%.1f\n", max_h)
end 


%%
%%QUESTION 5
%PROBABILITY
%%
clear all;
for k=1:5
    fprintf('\n');
    fprintf('\nQ5V%d\n', k);
    %question parameters
    lamb=0.01*randi([100, 250]);
    nstock=randi([3, 6]);
    %
    lamb5=lamb*5; %poisson mean for next 5 days
    %get the probability list.
    prob=zeros(1, nstock+1);
    for num=0:nstock-1
        temp=(exp(-lamb5)*(lamb5^num))/(factorial(num));
        prob(num+1)=temp;
    end
    temp=sum(prob);
    prob(nstock+1)=1-temp;
    exp=[];
    for rep=1:nstock+1
        exp(rep)=prob(rep)*(rep-1);
    end
    final_exp=sum(exp);
    %display the question
    fprintf("An island has a vending machine which gets repaired frequently. A repair tool can be used to repir it and a single tool can be used to repair it only once\n")
    fprintf("The number of such failures per day can be modelled as a poisson process with a mean of %0.3d \n", lamb)
    fprintf("A ship comes once every 5 days to restock %d repair tools. If all repair tools are utilised before the next ship arrives, the vending machine cant be used. \n", nstock)
    fprintf("Find the expected number of repairs which needs to be done on the vending machine ")
    %display the options:
    format short
    ops=[final_exp final_exp*0.8 final_exp*1.1 final_exp*0.5];
    choices=ops;
    choices_rand=choices(randperm(4));
    for i=1:4
        fprintf("\n %c.", i-1+'a');
        fprintf("%.2f", choices_rand(i));
        %disp(choices_rand(i));
        if(choices_rand(i)==choices(1))
            Ansopt=char(i-1+'a');
        end
    end

    fprintf('\n correct option is: option %c \n', Ansopt);
    %solution
    fprintf("\n")
    fprintf("SOLUTION  \n")
    format bank
    fprintf("The poisson parameter for 5 days cvan be calculated as 5 * %.2f =%.2f \n", lamb, lamb5)
    fprintf("the probability for each number of repair(n) is found to be exp(\x3BB)*\x3BB^(n)/(n!)  \n")
    fprintf("The probability of total number of repairs =total number of restocks can be found using the property \x2211 P(X)=1  \n")
    fprintf("E(X)=sum(P(X)*X)   \n")
    fprintf("multiplying the probabilities with corresponding event, we get the expected number =%.2f   \n", final_exp)
    clear
end
%%
%% QUESTION 6
%% STATISTICS
clear all
for i = 1:5
    %title
    fprintf('\n');
    fprintf('\nQ6V%d\n',i);
    
    %function
    temp=randi([1,2]);
    if temp==1
        %function = 1/(a+1-y)
        a = 0.01*randi(400, 1, 1);
        weights=randi([1, 10],1, 2);
        %correct answer
        sol=(1+a/2); %explanation is given in solution
        sol=weights(1)*sol +weights(2);
        %display the question (variant 1)
        fprintf(" \nY is a continuous uniform Random Variable on [0, a] where a=%.2f, Y~Unif[0, a])", a)
        fprintf("\n X=x follows an exponential random variable with mean 1/k where k=(a+1-y)")
        fprintf("\n If E(X) is the expected value of X, find %.0f*E(X) + %.2f\n ", weights(1), weights(2))
        ops=[sol sol-1 sol+1 sol-2];
        choices=ops;
        choices_rand=choices(randperm(4));
        for i=1:4
            fprintf("\n %c.", i-1+'a');
            fprintf("%.2f", choices_rand(i));
            %disp(choices_rand (i));
            if(choices_rand(i)==choices(1))
                Ansopt=char (i-1+'a');
            end
        end
        %correct option
        fprintf("\n")
        fprintf('\n correct option is: option %c \n', Ansopt);
        %SOLUTION
        fprintf("\nSOLUTION")
        fprintf("\nThe probability density function of Y is 1/a for a between 0, a and 0 elsewhere")
        fprintf("\nSince the mean is given to be 1/k, the parameter lambda=k")
        fprintf("\nThe expected value of X is calculated as integral 0 to a E[X|Y=y]f(y) where f(y) indicates probability density function of y")
        fprintf("\nThe quantity E(X|Y=y) =(1/lambda)=a+1-y. Substituting, we get E[X] as (1+ a/2)")
        fprintf("\nFor a=%.2f, Expected value= %.2f\n", a, 1+(a/2))
    else
        %function=1/(a^2 + y)
        a = 0.01*randi(400, 1, 1);
        weights=randi([1, 10],1, 2);
        %correct answer
        sol=(a^2 + a/2); %explanation given in solution
        sol=weights(1)*sol +weights(2);
        %display the question (variant 1)
        fprintf(" \nY is a continuous uniform Random Variable on [0, a],where a=%.2f Y~Unif[0, a])", a)
        fprintf("\nX=x follows an exponential random variable with mean k where k=(a^2 +y)")
        fprintf("\nIf E(X) is the expected value of X, find %.0f*E(X)+ %.2f\n ", weights(1), weights(2))
        ops=[sol sol-1 sol+1 sol-2];
        choices=ops;
        choices_rand=choices(randperm(4));
        for i=1:4
            fprintf("\n %c.", i-1+'a');
            fprintf("%.2f", choices_rand(i));
            %disp(choices_rand (i));
            if(choices_rand(i)==choices(1))
                Ansopt=char (i-1+'a');
            end
        end
        %correct option
        fprintf("\n")
        fprintf('\n correct option is: option %c \n', Ansopt);
        %SOLUTION
        fprintf("\nSOLUTION")
        fprintf("\nThe probability density function of Y is 1/a for a between 0, a and 0 elsewhere")
        fprintf("\nSince the mean is given to be k, the parameter lambda=1/k")
        fprintf("\nThe expected value of X is calculated as integral 0 to a E[X|Y=y]f(y) where f(y) indicates probability density function of y")
        fprintf("\nThe quantity E(X|Y=y) =(1/lambda)=a^2 + y. Substituting, we get E[X] as ((a^2)+ a/2)")
        fprintf("\nFor a=%.2f, Expected value= %.2f\n", a, a^2 +a/2)
    end
end 
%%
function [c, ceq] = nonlcon1(x)
    c = [];
    ceq = [];
end

%%
function [c, ceq] = nonlcon(x)
    c = [x(1)^2 + x(2)^2 - 1; 2*x(1) - x(2)^2];
    ceq = [];
end