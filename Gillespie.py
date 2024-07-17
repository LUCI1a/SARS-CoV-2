import numpy as np
import matplotlib.pyplot as plt

# 定义反应速率常数
k_bind,d_V,k_diss,k_fuse,k_uncoat,d_endosome,k_transl,f_ORF1,d_NSP,k_tr_neg,K_NSP,d_gRNA,d_gRNA_neg,k_tr_pos,k_complex,K_N,f_N,f_SP,d_N,d_SP,n_SP,n_N,K_Vrel,k_assemb,d_N_gRNA,k_release,d_assembled=12,0.12,0.61,0.5,0.5,0.06,45360,21000,0.069,3,100,0.2,0.1,1000,0.4,5000000,1200,10000,0.023,0.044,2000,456,1000,1,0.2,8,0.06

# 初始化系统状态

V_free,V_bound,V_endosome,gRNA_pos,NSP,gRNA_neg,gRNA,N,SP,N_gRNA,V_assembled,V_released = 100,0,0,0,0,0,0,0,0,0,0,0
total_time = 0
time_points = [total_time]
V_frees=[V_free]
V_bounds=[V_bound]
V_endosomes=[V_endosome]
gRNA_poss=[gRNA_pos]
NSPs=[NSP]
gRNA_negs=[gRNA_neg]
gRNAs=[gRNA]
Ns=[N]
SPs=[SP]
N_gRNAs=[N_gRNA]
V_assembleds=[V_assembled]
V_releaseds=[V_released]
# 设置模拟时间
simulation_time = 5

def reaction_rate(V_free,V_bound,V_endosome,gRNA_pos,NSP,gRNA_neg,gRNA,N,SP,N_gRNA,V_assembled,V_released):
    theta_RdRp = NSP / (NSP + K_NSP)
    theta_complex = N / (N + K_N)
    theta_assemb = SP / (SP + K_Vrel * n_SP)
    reaction_rate = np.zeros(23)
    reaction_rate[0] = V_free*k_bind
    reaction_rate[1] = V_free*d_N
    reaction_rate[2] = V_bound*k_diss
    reaction_rate[3] = V_bound*d_V
    reaction_rate[4] = V_bound*k_fuse
    reaction_rate[5] = V_endosome*k_uncoat
    reaction_rate[6] = V_endosome*d_endosome
    reaction_rate[7]= gRNA_pos*k_transl/f_ORF1
    reaction_rate[8] = NSP*d_NSP
    reaction_rate[9] = k_tr_neg*gRNA_pos*theta_RdRp
    reaction_rate[10] = gRNA_neg*d_gRNA_neg
    reaction_rate[11]= k_tr_pos*gRNA_neg*theta_RdRp
    reaction_rate[12]= gRNA*d_gRNA
    reaction_rate[13]= gRNA*k_complex*theta_complex
    reaction_rate[14]= gRNA*k_transl/f_N
    reaction_rate[15]= k_complex*n_N*theta_complex*gRNA+d_N*N
    reaction_rate[16]= gRNA*k_transl/f_SP
    reaction_rate[17]= k_assemb*n_SP*theta_assemb*N_gRNA+d_SP*SP
    reaction_rate[18]= k_assemb*theta_assemb*N_gRNA
    reaction_rate[19]= N_gRNA*d_N_gRNA
    reaction_rate[20]= V_assembled*k_release
    reaction_rate[21]= V_assembled*d_assembled
    reaction_rate[22]= V_released*d_V

    return (reaction_rate)
print(reaction_rate(V_free,V_bound,V_endosome,gRNA_pos,NSP,gRNA_neg,gRNA,N,SP,N_gRNA,V_assembled,V_released))
# Gillespie算法模拟
while total_time < simulation_time:
    # 计算反应速率
    reactionrate=reaction_rate(V_free,V_bound,V_endosome,gRNA_pos,NSP,gRNA_neg,gRNA,N,SP,N_gRNA,V_assembled,V_released)

    total_rate=np.sum(reactionrate[:])
    step1 = np.sum(reactionrate[:1]) / total_rate
    step2 = np.sum(reactionrate[:2]) / total_rate
    step3 = np.sum(reactionrate[:3]) / total_rate
    step4 = np.sum(reactionrate[:4]) / total_rate
    step5 = np.sum(reactionrate[:5]) / total_rate
    step6 = np.sum(reactionrate[:6]) / total_rate
    step7 = np.sum(reactionrate[:7]) / total_rate
    step8 = np.sum(reactionrate[:8]) / total_rate
    step9 = np.sum(reactionrate[:9]) / total_rate
    step10 = np.sum(reactionrate[:10]) / total_rate
    step11= np.sum(reactionrate[:11]) / total_rate
    step12= np.sum(reactionrate[:12]) / total_rate
    step13= np.sum(reactionrate[:13]) / total_rate
    step14= np.sum(reactionrate[:14]) / total_rate
    step15= np.sum(reactionrate[:15]) / total_rate
    step16= np.sum(reactionrate[:16]) / total_rate
    step17= np.sum(reactionrate[:17]) / total_rate
    step18= np.sum(reactionrate[:18]) / total_rate
    step19= np.sum(reactionrate[:19]) / total_rate
    step20= np.sum(reactionrate[:20]) / total_rate
    step21= np.sum(reactionrate[:21]) / total_rate
    step22= np.sum(reactionrate[:22]) / total_rate

    rand1 = np.random.rand()
    rand2 = np.random.rand()

    # 计算时间步长
    delta_t = -np.log(rand1) / total_rate

    # 选择发生的反应
    if rand2 < step1:
        V_free -= 1
        V_bound += 1
    elif rand2>=step1 and rand2<step2:
        V_free -= 1
    elif rand2>=step2 and rand2<step3:
        V_bound -= 1
        V_free += 1
    elif rand2>=step3 and rand2<step4:
        V_bound -= 1
    elif rand2>=step4 and rand2<step5:
        V_bound -= 1
        V_endosome += 1
    elif rand2>=step5 and rand2<step6:
        V_endosome -= 1
        gRNA_pos += 1
    elif rand2>=step6 and rand2<step7:
        V_endosome -= 1
    elif rand2>=step7 and rand2<step8:
        NSP += 1
        #gRNA_pos -= 1
    elif rand2>=step8 and rand2<step9:
        NSP -= 1
    elif rand2>=step9 and rand2<step10:
        gRNA_neg += 1
    elif rand2>=step10 and rand2<step11:
        gRNA_neg -= 1
    elif rand2 >= step11 and rand2 < step12:
        gRNA += 1
        #gRNA_neg -= 1
    elif rand2 >= step12 and rand2 < step13:
        gRNA -= 1
    elif rand2 >= step13 and rand2 < step14:
        gRNA -= 1
        N_gRNA += 1
    elif rand2 >= step14 and rand2 < step15:
        N += 1
    elif rand2 >= step15 and rand2 < step16:
        N -= 1
    elif rand2 >= step16 and rand2 < step17:
        SP += 1
    elif rand2 >= step17 and rand2 < step18:
        SP -= 1
    elif rand2 >= step18 and rand2 < step19:
        N_gRNA -= 1
        V_assembled += 1
    elif rand2 >= step19 and rand2 < step20:
        N_gRNA -= 1
    elif rand2 >= step20 and rand2 < step21:
        V_assembled -= 1
        V_released += 1
    elif rand2 >= step21 and rand2 < step22:
        V_assembled -= 1
    else:
        V_released -= 1
# 更新总时间和系统状态
    total_time += delta_t
    time_points.append(total_time)
    print(total_time)
    V_frees.append(V_free)
    V_bounds.append(V_bound)
    V_endosomes.append(V_endosome)
    gRNA_poss.append(gRNA_pos)
    NSPs.append(NSP)
    gRNA_negs.append(gRNA_neg)
    gRNAs.append(gRNA)
    Ns.append(N)
    SPs.append(SP)
    N_gRNAs.append(N_gRNA)
    V_assembleds.append(V_assembled)
    V_releaseds.append(V_released)
  # 绘制模拟结果
# plt.plot(time_points, V_frees, label='Vfree')
# plt.plot(time_points, V_bounds, label='Vbound')
# plt.plot(time_points, V_endosomes, label='Vendosome')
# plt.plot(time_points, gRNA_poss, label='gRNA(+)')
# plt.plot(time_points, NSPs, label='NSP')
# plt.plot(time_points, gRNA_negs, label='gRNA(-)')
# plt.plot(time_points, gRNAs, label='gRNA')
# plt.plot(time_points, Ns, label='N')
# plt.plot(time_points, SPs, label='SP')
# plt.plot(time_points, N_gRNAs, label='N-gRNA')
plt.plot(time_points, V_assembleds, label='Vassembled')
plt.plot(time_points, V_releaseds, label='Vreleased')

plt.xlabel('Time')
plt.ylabel('Molecule Count')
plt.legend()
plt.show()
