function [value, isterminal, direction] = add_vemurafenib(t, y)

t_add = 5e4;

value = [sign(t - t_add)];
isterminal = [1];
direction = [+1];

end