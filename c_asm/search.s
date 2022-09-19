%define cost_best       8
%define cost_concat_1   16
%define cost_concat_2   24
%define rev_seq_cost    32
%define cost_new        40
%define solut           48
%define info            56
%define I               64
%define J               72

struc   tSeqInfo
    .T: resq    1
    .C: resq    1
    .W: resq    1
endstruc

struc  tSolution
    .seq:       resq    1
    .cost:      resq    1
    .s:         resq    1
    .s_size:    resd    1
endstruc

struc tInfo
    .cost:      resq    1
    .dimen:     resd    1
endstruc

section .text

    extern subseq_load
    global search_two_opt_asm

search_two_opt_asm:

    push    r15
    push    r14
    push    r13
    push    r12
    push    rbx
    push    rbp
    mov     rbp, rsp

    xor     rbx, rbx

    push    rbx     ; cost_best
    push    rbx     ; concat_1
    push    rbx     ; concat_2
    push    rbx     ; rev_seq_cost
    push    rbx     ; cost_new
    push    rdi     ; tSolution * solut
    push    rsi     ; tInfo * solut
    push    rbx     ; I
    push    rbx     ; J

    ; rbx -> i
    ; rcx -> j

    ; rdx -> i_prev
    ; r10 -> j_next

    ; r14 -> i_bound
    ; r15 -> j_bound

    ; xmm0 -> cost_new
    ; xmm1 -> cost_concat_1
    ; xmm2 -> cost_concat_2
    ; xmm3 -> rev_seq_cost


    xorpd     xmm3, xmm3
    movsd   xmm7, QWORD[inf]

    mov     rbx, 1                          ; i = 1
    mov     r14, QWORD [rbp - info]
    mov     r14d, DWORD [r14 + tInfo.dimen]
    dec     r14d                             ; i_bound = dimen-1

    jmp     for_i_cond
for_i:
    mov     rdx, rbx
    dec     rdx                             ; i_prev = i-1

    ; rev_seq_cost = seq[i][i+1].T (xmm3)
    mov     r11, rdi
    ;mov     r11, QWORD [rbp - solut]
    mov     r11, QWORD [r11 + tSolution.seq] ; r11 -> seq
    mov     r11, QWORD [r11 + rbx*8]         ; r11 -> seq[i]
    mov     r12, rbx
    inc     r12
    imul    r12, 24  ; r12 -> offset
    lea     r11, [r11 + r12]
    movsd   xmm3, QWORD [r11 + tSeqInfo.T]



    ; FOR J
    mov     rcx, rbx        ; j = i+2
    add     rcx, 2

    mov     r15d,r14d       ; j_bound = dimen
    inc     r15d
    jmp for_j_cond
for_j:
    mov     r10, rcx        ; j_next = j+1
    inc     r10 

    mov     r11, QWORD [rdi + tSolution.s]
    mov     r12d, DWORD [r11 + rcx*4 - 4]          ; s[j-1]

    mov     r8, QWORD [rsi + tInfo.cost]
    mov     r8, QWORD [r8 + r12*8]              ; cost[s[j-1]]
    mov     r12d, DWORD [r11 + rcx*4]              ; s[j]
    movsd     xmm4, QWORD [r8 + r12*8]              ; cost[s[j-1]][s[j]]
    ; xmm4 -> info->cost[solut->s[j-1]][solut->s[j]]

    mov     r11, QWORD [rdi + tSolution.seq]
    mov     r11, QWORD [r11 + rbx*8]            ; solut->seq[i]
    mov     r8, rcx
    imul     r8, 24
    lea     r11, [r11 + r8]            ; solut->seq[i][j]
    movsd   xmm5, QWORD [r11 + tSeqInfo.W]
    subsd   xmm5, [um]            ; (solut->seq[i][j].W - 1.0)
     
    mulsd   xmm4, xmm5
    addsd   xmm3, xmm4
    
    ; cost_concat_1 
    mov     r11, QWORD [rdi + tSolution.seq]
    mov     r11, QWORD [r11]            ; solut->seq[0]
    mov     r8, rdx
    imul    r8, 24
    lea     r11, [r11 + r8]            ; solut->seq[0][i_prev]
    movsd   xmm1, QWORD [r11 + tSeqInfo.T]  ; cost_concat_1 = solut->seq[0][i_prev].T
    
    mov     r11,  QWORD [rdi + tSolution.s]
    mov     r12d, DWORD [r11 + rcx*4]           ; s[j]
    mov     r8,   QWORD [rsi + tInfo.cost]
    mov     r8,   QWORD [r8 + r12*8]              ; cost[s[j]]
    mov     r12d, DWORD [r11 + rdx*4]              ; s[i_prev]
    movsd   xmm4, QWORD [r8 + r12*8]              ; cost[s[j]][s[i_prev]]
    addsd   xmm1, xmm4                            ; cost_concat_1 += cost[s[j]][s[i_prev]]
    

    ; cost_concat_2 
    mov     r11, QWORD [rdi + tSolution.seq]
    mov     r11, QWORD [r11 + rbx*8]            ; solut->seq[i]
    mov     r8, rcx
    imul    r8, 24
    lea     r11, [r11 + r8]            ; solut->seq[i][j]
    movsd   xmm2, QWORD [r11 + tSeqInfo.T]  ; cost_concat_2 = solut->seq[i][j].T
    
    mov     r11,  QWORD [rdi + tSolution.s]
    mov     r12d, DWORD [r11 + r10*4]           ; s[j_next]
    mov     r8,   QWORD [rsi + tInfo.cost]
    mov     r8,   QWORD [r8 + r12*8]              ; cost[s[j_next]]
    mov     r12d, DWORD [r11 + rbx*4]              ; s[i]
    movsd   xmm4, QWORD [r8 + r12*8]              ; cost[s[j_next]][s[i]]
    addsd   xmm2, xmm4                            ; cost_concat_2 += cost[s[j_next]][s[i]]
    addsd   xmm2, xmm1                            ; cost_concat_2 += cost_concat_1
    

    ; solut->seq[i][j].W * cost_concat_1
    mov     r11, QWORD [rdi + tSolution.seq]
    mov     r11, QWORD [r11 + rbx*8]            ; solut->seq[i]
    mov     r8, rcx
    imul    r8, 24
    lea     r11, [r11 + r8]                     ; solut->seq[i][j]
    movsd   xmm4, QWORD [r11 + tSeqInfo.W]      ; solut->seq[i][j].W
    mulsd   xmm1, xmm4


    ; solut->seq[j_next][info->dimen].W * cost_concat_2
    mov     r11, QWORD [rdi + tSolution.seq]
    mov     r11, QWORD [r11 + r10*8]            ; solut->seq[j_next]
    mov     r8d, DWORD [rsi + tInfo.dimen]
    imul    r8, 24
    lea     r11, [r11 + r8]                     ; solut->seq[j_next][info->dimen]
    movsd   xmm4, QWORD [r11 + tSeqInfo.W]      ; solut->seq[j_next][info->dimen].W
    movsd   xmm5, QWORD [r11 + tSeqInfo.C]      ; solut->seq[j_next][info->dimen].C
    mulsd   xmm2, xmm4

    ; solut->seq[0][i_prev].C
    mov     r11, QWORD [rdi + tSolution.seq]
    mov     r11, QWORD [r11]            ; solut->seq[0]
    mov     r8, rdx
    imul    r8, 24
    lea     r11, [r11 + r8]                     ; solut->seq[0][i_prev]
    movsd   xmm4, QWORD [r11 + tSeqInfo.C]      ; solut->seq[0][i_prev].C


    movsd   xmm0, xmm4                          ; cost_new = solut->seq[0][i_prev].C
    addsd   xmm0, xmm5                          ; cost_new += solut->seq[j_next][info->dimen].C
    addsd   xmm0, xmm1                          ; cost_new += solut->seq[i][j].W * cost_concat_1
    addsd   xmm0, xmm2                          ; cost_new += solut->seq[j_next][info->dimen].W * cost_concat_2
    addsd   xmm0, xmm3                          ; cost_new += rev_seq_cost

    

    movsd   xmm4, xmm0
    comisd   xmm7, xmm4       ; cmp lt
    jbe      for_j_if_false      ; if (cost_new < cost_best)
    movsd   xmm7, xmm0          ; cost_best = cost_new
    mov     DWORD [rbp-I], ebx
    mov     DWORD [rbp-J], ecx
    for_j_if_false:



    inc     rcx
for_j_cond:
    cmp     ecx, r15d
    jl      for_j
    

    ; FOR_I conditions
    inc     rbx                             ; i++
for_i_cond:
    cmp     rbx, r14
    jl      for_i

    ; if (cost_best < solut->cost)
    movsd     xmm4, QWORD [rdi + tSolution.cost]
    comisd  xmm4, xmm7
    jbe two_opt_if_false

    mov     eax, DWORD [rbp-I] 
    add     eax, DWORD [rbp-J] 
    mov     bx,  2
    cdq
    div     bx

    mov     ebx, DWORD [rbp-I]
    mov     ecx, DWORD [rbp-J]
    mov     edx, eax

    mov     r11, QWORD [rdi+tSolution.s]

    jmp     for_rev_seq_cond
for_rev_seq:

    mov     eax, DWORD [r11 + rbx*4]
    mov     r8d, DWORD [r11+rcx*4]
    mov     DWORD [r11+rbx*4], r8d
    mov     DWORD [r11+rcx*4], eax


    inc     rbx
    dec     rcx
for_rev_seq_cond:
    cmp     ebx, edx
    jle     for_rev_seq


    mov rax, 0
    call subseq_load wrt ..plt
    mov rax, 1
    jmp two_opt_return

two_opt_if_false:
    mov rax, 0

two_opt_return:

    pop     rbx
    pop     rbx
    pop     rbx
    pop     rbx
    pop     rbx
    pop     rbx
    pop     rbx
    pop     rbx
    pop     rbx

    pop    rbp
    pop    rbx
    pop    r12
    pop    r13
    pop    r14
    pop    r15

    ret

section .data
    um: dq 0x3ff0000000000000
    inf: dq 0x7FF0000000000000 
