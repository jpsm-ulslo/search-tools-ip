import os
import time
import shutil
import pandas as pd
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import TimeoutException
from selenium.webdriver.common.action_chains import ActionChains

from _utils import iniciar_chrome, digitar_simulado, aguardar_download

DOWNLOAD_DIR = os.getcwd()

def distribuidores_ativos(caminho_xlsx, verbose=False):
    """
    Lê ficheiro Excel do Infarmed e devolve lista única de distribuidores
    com pelo menos um CDM atualmente comercializado.
    
    Args:
        caminho_xlsx (str): Caminho para o Excel exportado pelo Infarmed
        verbose (bool): Se True, imprime informação de debug
    
    Returns:
        list[str]: distribuidores únicos ordenados alfabeticamente
    """
    # Carregar Excel
    xls = pd.ExcelFile(caminho_xlsx)
    df = pd.read_excel(caminho_xlsx, sheet_name=xls.sheet_names[0]).fillna("")

    if verbose:
        print("Sheets encontradas:", xls.sheet_names)
        print("Colunas encontradas:", df.columns.tolist())

    # Normalizar nomes de colunas
    colunas_normalizadas = {c.strip().lower(): c for c in df.columns}

    # Identificar colunas
    col_estado_cdm = next((col for key, col in colunas_normalizadas.items() if "estado" in key and "distribuidor" not in key), None)
    col_estado_dist = next((col for key, col in colunas_normalizadas.items() if "estado" in key and "distribuidor" in key), None)
    col_distribuidor = colunas_normalizadas.get("distribuidor")

    if not col_estado_cdm or not col_estado_dist or not col_distribuidor:
        if verbose:
            print("⚠️ Não foram encontradas todas as colunas necessárias.")
        return []

    # Filtrar apenas CDMs e distribuidores ativos
    df_filtrado = df[
        df[col_estado_cdm].str.strip().str.lower().str.contains("comercializado")
        & df[col_estado_dist].str.strip().str.lower().str.contains("comercializado")
    ]

    # Extrair lista única de distribuidores
    distribuidores_unicos = sorted(df_filtrado[col_distribuidor].dropna().unique())

    if verbose:
        print(f"{len(distribuidores_unicos)} distribuidores únicos encontrados.")

    return distribuidores_unicos

def ler_cdms_excel(caminho_xlsx):
    """Ler arquivo Excel exportado pelo Infarmed e retornar lista de CDMs com Estado = Comercializado."""
    cdms = []

    # Verificar sheets disponíveis
    xls = pd.ExcelFile(caminho_xlsx)
    print("Sheets encontradas:", xls.sheet_names)

    # Tentar abrir a primeira sheet (normalmente é a que tem os resultados)
    df = pd.read_excel(caminho_xlsx, sheet_name=xls.sheet_names[0])
    print("Colunas encontradas:", df.columns.tolist())

    # Normalizar nomes das colunas
    colunas_normalizadas = {c.strip().lower(): c for c in df.columns}

    # Procurar nome da coluna que contenha "estado"
    col_estado = next((col for key, col in colunas_normalizadas.items() if "estado" in key), None)
    if not col_estado:
        print("⚠️ Coluna de Estado não encontrada.")
        return cdms

    # Iterar sobre linhas
    for _, linha in df.iterrows():
        estado = str(linha.get(col_estado, "")).strip().lower()
        if "comercializado" in estado:
            cdms.append({
                "Referência": str(linha.get(colunas_normalizadas.get("referência", "Referência"), "")).strip(),
                "Fabricante": str(linha.get(colunas_normalizadas.get("fabricante", "Fabricante"), "")).strip(),
                "CDM": str(linha.get(colunas_normalizadas.get("cdm", "CDM"), "")).strip(),
                "Distribuidor": str(linha.get(colunas_normalizadas.get("distribuidor", "Distribuidor"), "")).strip(),
            })

    print(f"{len(cdms)} CDMs Comercializados encontrados.")
    return cdms


def aguardar_e_renomear_download(download_dir, cdm, timeout=60):
    """
    Aguarda pelo download terminar e renomeia para {CDM}_folheto.pdf.
    Detecta o arquivo novo na pasta de downloads.
    """
    fim = time.time() + timeout
    arquivos_antes = set(os.listdir(download_dir))

    while time.time() < fim:
        arquivos_atuais = set(os.listdir(download_dir))
        novos_arquivos = arquivos_atuais - arquivos_antes
        for nome in novos_arquivos:
            if nome.endswith(".crdownload"):
                continue  # Download ainda em andamento
            if nome.endswith(".pdf"):
                antigo_caminho = os.path.join(download_dir, nome)
                novo_caminho = os.path.join(download_dir, f"{cdm}_folheto.pdf")
                shutil.move(antigo_caminho, novo_caminho)
                return True
        time.sleep(0.5)
    return False


def processar_cdms_individualmente(cdms_info, relatorio_path="relatorio_cdms.xlsx"):
    """
    Para cada CDM:
    - Preencher pesquisa, abrir modal, extrair dados
    - Baixar folheto informativo (com validação)
    - Guardar relatório final em Excel/CSV
    """
    driver = iniciar_chrome()
    wait = WebDriverWait(driver, 30)
    resultados = []

    try:
        driver.get("https://www.infarmed.pt/web/infarmed/pesquisa-dispositivos")

        for idx, info in enumerate(cdms_info, 1):
            cdm = info["CDM"]
            resultado = {"CDM": cdm, "Status": "N/A", "Folheto_Baixado": False}
            print(f"[{idx}/{len(cdms_info)}] Processando CDM: {cdm}")

            try:
                # Campo de pesquisa
                cdm_input = wait.until(EC.presence_of_element_located(
                    (By.ID, "_SIDMPesquisaDispositivos_WAR_SIDMportlet_:pesquisaForm:infarmed5_input")
                ))
                cdm_input.clear()
                digitar_simulado(cdm_input, cdm)

                # Botão pesquisar
                search_btn = wait.until(EC.element_to_be_clickable(
                    (By.ID, "_SIDMPesquisaDispositivos_WAR_SIDMportlet_:pesquisaForm:pesquisarBtn")
                ))
                driver.execute_script("arguments[0].click();", search_btn)
                print("  - Pesquisa realizada")

                # Ícone de detalhes
                detalhes_icon = wait.until(EC.element_to_be_clickable(
                    (By.CSS_SELECTOR, "div.icon-file-alt[title='Detalhes']")
                ))
                driver.execute_script("arguments[0].click();", detalhes_icon)
                print("  - Modal de detalhes aberto")

                # Modal
                modal = wait.until(EC.visibility_of_element_located(
                    (By.ID, "_SIDMPesquisaDispositivos_WAR_SIDMportlet_:detalheForm:detalheDialog")
                ))

                # Extrair alguns campos principais
                campos = [
                    "Referência", "Designação/Nome Comercial", "Tipo de Dispositivo",
                    "Marca", "Modelo", "Classe", "Família", "Estado de Comercialização",
                    "Código NPDM", "Fabricante"
                ]
                for campo in campos:
                    try:
                        elemento = modal.find_element(By.XPATH, f".//span[text()='{campo}']/following-sibling::*[1]")
                        resultado[campo] = elemento.text.strip()
                    except:
                        resultado[campo] = None

                # Baixar folheto com salvaguarda
                resultado["Folheto_Baixado"] = baixar_folheto(driver, modal, cdm)

                resultado["Status"] = "Processado"

                # Fechar modal
                try:
                    botao_fechar = wait.until(EC.element_to_be_clickable(
                        (By.CSS_SELECTOR, "a.ui-dialog-titlebar-close")
                    ))
                    driver.execute_script("arguments[0].click();", botao_fechar)
                    wait.until(EC.invisibility_of_element_located(
                        (By.ID, "_SIDMPesquisaDispositivos_WAR_SIDMportlet_:detalheForm:detalheDialog")
                    ))
                    print(f"  - Modal fechado para CDM {cdm}")
                except TimeoutException:
                    print(f"  - Não foi possível fechar o modal para CDM {cdm}")

            except TimeoutException:
                resultado["Status"] = "Não encontrado"
                print(f"  - CDM {cdm} não encontrado ou sem resultados.")

            resultados.append(resultado)

    finally:
        driver.quit()
        print("Chrome fechado após processar CDMs.")

        # Guardar relatório final em Excel
        df = pd.DataFrame(resultados)
        df.to_excel(relatorio_path, index=False)
        print(f"📄 Relatório salvo em: {relatorio_path}")

def baixar_folheto(driver, modal, cdm, max_tentativas=3):
    """
    Baixa o folheto informativo de um CDM com salvaguardas:
    - Espera o download terminar antes de prosseguir
    - Retenta caso falhe
    """
    for tentativa in range(1, max_tentativas + 1):
        try:
            folheto_icon = modal.find_element(By.CSS_SELECTOR,
                "span.icon-download-alt[title='Folheto Informativo']")
            driver.execute_script("arguments[0].click();", folheto_icon)
            print(f"  - Tentativa {tentativa}: download iniciado para CDM {cdm}")

            if aguardar_e_renomear_download(DOWNLOAD_DIR, cdm, timeout=60):
                caminho_pdf = os.path.join(DOWNLOAD_DIR, f"{cdm}_folheto.pdf")
                if os.path.exists(caminho_pdf) and os.path.getsize(caminho_pdf) > 0:
                    print(f"  - PDF salvo e validado: {caminho_pdf}")
                    return True
                else:
                    print(f"  - Arquivo inválido ou vazio: {caminho_pdf}")
            else:
                print(f"  - Timeout: PDF do CDM {cdm} não encontrado")

        except Exception as e:
            print(f"  - Erro ao tentar baixar PDF do CDM {cdm}: {e}")

        print("  - Retentando download...")

    print(f"❌ Falha final: não foi possível obter folheto do CDM {cdm}")
    return False
