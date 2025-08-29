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


def processar_cdms_individualmente(cdms_info):
    """
    Para cada CDM:
    - Preencher pesquisa, abrir modal, extrair dados
    - Baixar folheto informativo, esperar download terminar
    - Fechar modal
    """
    driver = iniciar_chrome()
    wait = WebDriverWait(driver, 30)

    try:
        driver.get("https://www.infarmed.pt/web/infarmed/pesquisa-dispositivos")

        for idx, info in enumerate(cdms_info, 1):
            cdm = info["CDM"]
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

                # Extrair campos
                campos = [
                    "CDM", "Referência", "Designação/Nome Comercial", "Tipo de Dispositivo",
                    "Marca", "Modelo", "Classe", "Família", "Estado de Comercialização",
                    "Código NPDM", "Fabricante"
                ]
                for campo in campos:
                    try:
                        elemento = modal.find_element(By.XPATH, f".//span[text()='{campo}']/following-sibling::*[1]")
                        info[campo] = elemento.text.strip()
                    except:
                        info[campo] = None

                # Baixar folheto informativo
                info["Folheto_Baixado"] = False
                try:
                    folheto_icon = modal.find_element(By.CSS_SELECTOR,
                        "span.icon-download-alt[title='Folheto Informativo']")
                    driver.execute_script("arguments[0].click();", folheto_icon)
                    print(f"  - Folheto iniciado para CDM {cdm}")

                    # Esperar download terminar e renomear
                    if aguardar_e_renomear_download(DOWNLOAD_DIR, cdm):
                        info["Folheto_Baixado"] = True
                        print(f"  - PDF renomeado para {cdm}_folheto.pdf")
                    else:
                        print(f"  - Timeout: PDF do CDM {cdm} não encontrado")
                except:
                    print(f"  - Folheto não encontrado para CDM {cdm}")

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
                print(f"  - CDM {cdm} não encontrado ou sem resultados.")

    finally:
        driver.quit()
        print("Chrome fechado após processar CDMs.")
