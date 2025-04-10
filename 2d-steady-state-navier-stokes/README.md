Código para escoamento laminar em canal em regime permanente, com $\text{Re} = 10.0$.

$$
\text{Re} \mathbf{u} \cdot \nabla \mathbf{u} = - \nabla p + \nabla^2 \mathbf{u}
$$

Condições de contorno (explícitas):
- Velocidade uniforme na entrada com $\mathbf{u}=(10.0,0)$
- Velocidade nula nas paredes
- Pressão igual a zero na saída

A componente vertical da velocidade tem comportamento anômalo no início do domínio, fazendo com que as linhas de corrente vão em direção às paredes em vez de irem em direção ao centro do escoamento. Esse comportamento parece ter relação com a resolução da malha, fazendo-se presente e bem pronunciada quando a malha possui $n \times n$ elementos e desaparecendo quando é mais refinada na direção do escoamento (aqui, alinhado ao eixo $x$).

![](/2d-steady-state-navier-stokes/speed-profile-pictures/speed-profile-50x50.png)
![](/2d-steady-state-navier-stokes/speed-profile-pictures/speed-profile-100x100.png)
![](/2d-steady-state-navier-stokes/speed-profile-pictures/speed-profile-250x250.png)
![](/2d-steady-state-navier-stokes/speed-profile-pictures/speed-profile-250x50.png)

Para que a simulação com uma malha $500 \times 100$ convirja, é necessário aumentar o grau dos polinômios interpoladores de 2 para 3. Não há muito ganho de precisão quando comparado à malha de $250 \times 50$, mas o perfil da componente vertical fica mais suave.

![](/2d-steady-state-navier-stokes/speed-profile-pictures/speed-profile-500x100.png)

.