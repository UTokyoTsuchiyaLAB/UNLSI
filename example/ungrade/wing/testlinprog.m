% ランダムな係数行列とベクトル生成
rng(42); % 乱数の再現性を確保するためのシード
m = 5;    % 制約の数
n = 3;    % 変数の数
A = randn(m, n);
b = randn(m, 1);
c = randn(n, 1);

% linprogで解く
options_linprog = optimoptions('linprog', 'Display', 'off');
tic;
[x_linprog, fval_linprog, exitflag_linprog, output_linprog] = linprog(c, -A, -b, [], [], zeros(n,1), []);
time_linprog = toc;

% 自作の内点法で解く
tic;
[x_interior, lambda, fval_interior] = interiorPointInequalitySolverWithLagrange(A, b, c);
time_interior = toc;

% 結果の表示
disp('linprog result:');
disp(['Optimal x: ' num2str(x_linprog')]);
disp(['Optimal fval: ' num2str(fval_linprog)]);
disp(['Exit flag: ' num2str(exitflag_linprog)]);
disp(['Elapsed time: ' num2str(time_linprog) ' seconds']);

disp('Interior point solver result:');
disp(['Optimal x: ' num2str(x_interior')]);
disp(['Optimal fval: ' num2str(fval_interior)]);
disp(['Elapsed time: ' num2str(time_interior) ' seconds']);