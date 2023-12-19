# import sympy

# def get_formula(n):
#     formula_str = "1/4+2*(0"
#     for k in range(1,n):
#         formula_str += f"+8*{k}/(256+{k}**2)"
#     formula_str += ")+1/5"
#     return formula_str

# p2str = "1/16*("+get_formula(8)+")"
# print(p2str)
# print(sympy.simplify(p2str)) # wrong

ans = 0
for k in range(1,8):
    ans += 8*k/(256+k**2)
ans *= 2
sum1 = ans
ans += 1/4
ans += 1/5
ans *= 1/16
print(ans)

ans2 = 0
for k in range(8):
    ans2 += (32*k+16)/(4*k*k+4*k+1025)
ans2 *= 4
ans2 += 1/4+1/5 + sum1
ans2 *= 1/48
print(ans2)
