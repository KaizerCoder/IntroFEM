import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

def prepara_plot(func,var,x0: float,xf: float,label=None,pontos=100):
    """Prepara uma função para plotar, após todas as funções terem sido preparadas, é só colocar os labels dos eixos(opcional) e dar o comando plt.show() :D
    Args:
        func (_type_): Função a preparar (Sympy.symbols)
        var (_type_): Varivel dependente da função (Sympy.symbols)
        x0 (_type_): Ponto inicial da variavel dependente
        xf (_type_): Ponto Final da variavel dependente
        label (_type_, optional): _description_. Label da Função, se não preenchido, fico como a equação da função
        pontos (int, optional): _description_. Pontos, padrão 100.
    """    
    if label is None:
        label = sp.simplify(func)
    
    T_func = sp.lambdify(var, func, modules='numpy')

    # Gerar valores e plotar
    x_vals = np.linspace(x0, xf, pontos)
    y_vals = T_func(x_vals)

    plt.plot(x_vals, y_vals, label=f'{label}')
    
    
    
    
def verifica_condicoes(f, x, condicoes, verbose=False, tol=1e-6):
    """
    Verifica condições de valor e derivada em pontos específicos.

    Args:
        f: função simbólica f(x)
        x: variável simbólica
        condicoes: lista de tuplas do tipo:
            - ('valor', x0, valor_esperado)
            - ('derivada', x0, valor_esperado)
        verbose: exibe relatório no console
        tol: tolerância para considerar igual

    Returns:
        Lista de resultados [(ok: bool, erro_abs: float, mensagem: str)]
    """
    resultados = []

    for tipo, ponto, valor_esperado in condicoes:
        if tipo == 'valor':
            resultado = f.subs(x, ponto).evalf()
        elif tipo == 'derivada':
            resultado = sp.diff(f, x).subs(x, ponto).evalf()
        else:
            raise ValueError(f"Tipo desconhecido: {tipo}")
        
        erro = abs(resultado - valor_esperado)
        ok = erro < tol
        # msg = f"{tipo} em x={ponto}: esperado={valor_esperado}, obtido={resultado}, erro={erro:.2e}"
        # if verbose:
            # print(("✅" if ok else "❌") + " " + msg)
        # resultados.append((ok, erro, msg))
    
    return ok

def compatibiliza_cc_com_derivada(expr, x, variaveis, condicoes):
    """
    Ajusta a expressão para satisfazer uma ou mais condições (valores ou derivadas)
    eliminando variáveis livres.

    Args:
        expr: expressão simbólica (com coeficientes livres, ex: a0, a1, ...)
        x: variável simbólica independente
        variaveis: lista de variáveis livres (ex: [a0, a1, a2, a3])
        condicoes: lista de tuplas:
            - ('valor', x0, v0)
            - ('derivada', x1, v1)
            - ('2derivada', x1, v1)
            - ('3derivada', x1, v1)

    Returns:
        nova expressão com as variáveis eliminadas e dicionário de substituições
    """
    eqs = []

    for tipo, x0, valor in condicoes:
        if tipo == 'valor':
            eqs.append(sp.Eq(expr.subs(x, x0), valor))
        elif tipo == 'derivada':
            derivada = sp.diff(expr, x)
            eqs.append(sp.Eq(derivada.subs(x, x0), valor))
        elif tipo == '2derivada':
            derivada = sp.diff(expr, x,2)
            eqs.append(sp.Eq(derivada.subs(x, x0), valor))
        elif tipo == '3derivada':
            derivada = sp.diff(expr, x,3)
            eqs.append(sp.Eq(derivada.subs(x, x0), valor))
        else:
            raise ValueError(f"Tipo de condição não reconhecido: {tipo}")

    # Resolver as condições para um subconjunto das variáveis
    sol = sp.solve(eqs, variaveis, dict=True)

    if not sol:
        raise ValueError("Não foi possível compatibilizar as condições")

    # Substituir na expressão original
    expr_compat = expr.subs(sol[0])
    return expr_compat.simplify(),sol


def VetorColuna(vec):
    return sp.Matrix(vec).T.T
    
        