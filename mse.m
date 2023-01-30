function MSE = mse(target,estim_target)

if size(target,2)>1
    MSE = vec(target-estim_target)'*vec(target-estim_target);
else
    MSE = (target-estim_target)'*(target-estim_target);
end