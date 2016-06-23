function [  Violation_]=  burgersAdvection(dudt, lambda)

tvdFun = @(u) sum([abs(diff(u)); abs((u(1)-u(end)))]);

Violation_ = nan(1, length(lambda));
nSteps = 5;
for i = 1:numel(lambda);
    dt = lambda(i)*dudt.dfdx.dx;
    
    dudt.resetInitCondition();
    [~, y_] = dudt.getState();
    TV_ = nan(1, nSteps);
    TV_(1) = tvdFun(y_);
    
    for tt = 2:nSteps
        dudt.takeStep(dt);
        [~, y_] = dudt.getState();
        TV_(tt) = tvdFun(y_);
    end
    Violation_(i) = max([diff(TV_),1e-15]);
end

end