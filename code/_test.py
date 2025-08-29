import os
import time
import pandas as pd
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import TimeoutException

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


def processar_cdms_individualmente(cdms_info):
    """
    Para cada CDM, acessar a página de pesquisa do Infarmed, preencher o CDM,
    pesquisar e abrir diretamente o modal de detalhes, ignorando a tabela.
    """
    driver = iniciar_chrome()
    wait = WebDriverWait(driver, 20)

    try:
        driver.get("https://www.infarmed.pt/web/infarmed/pesquisa-dispositivos")

        for idx, info in enumerate(cdms_info, 1):
            cdm = info["CDM"]
            print(f"[{idx}/{len(cdms_info)}] Processando CDM: {cdm}")

            try:
                # 1. Campo de pesquisa
                cdm_input = wait.until(EC.presence_of_element_located(
                    (By.ID, "_SIDMPesquisaDispositivos_WAR_SIDMportlet_:pesquisaForm:infarmed5_input")
                ))
                cdm_input.clear()
                digitar_simulado(cdm_input, cdm)
                time.sleep(0.2)

                # 2. Botão de pesquisar
                search_btn = wait.until(EC.element_to_be_clickable(
                    (By.ID, "_SIDMPesquisaDispositivos_WAR_SIDMportlet_:pesquisaForm:pesquisarBtn")
                ))
                driver.execute_script("arguments[0].click();", search_btn)
                print("  - Pesquisa realizada")

                # 3. Esperar pelo ícone de detalhes do CDM (div.icon-file-alt) de forma confiável
                detalhes_icon = wait.until(EC.element_to_be_clickable(
                    (By.CSS_SELECTOR, "div.icon-file-alt[title='Detalhes']")
                ))
                driver.execute_script("arguments[0].scrollIntoView(true);", detalhes_icon)
                driver.execute_script("arguments[0].click();", detalhes_icon)
                print("  - Modal de detalhes aberto")

                # 4. Aguardar modal carregar
                modal = wait.until(EC.visibility_of_element_located(
                    (By.ID, "_SIDMPesquisaDispositivos_WAR_SIDMportlet_:detalheForm:detalheDialog")
                ))
                time.sleep(0.2)

                # 5. Extrair campos do modal
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

                print(f"  - Dados extraídos para CDM {cdm}: {info}")

            except TimeoutException:
                print(f"  - CDM {cdm} não encontrado ou modal não abriu.")

            # 6. Fechar modal antes de ir para próximo CDM
            try:
                close_btn = modal.find_element(By.CSS_SELECTOR, "a.ui-dialog-titlebar-icon.ui-dialog-titlebar-close")
                driver.execute_script("arguments[0].click();", close_btn)
                time.sleep(0.2)
            except:
                pass

    finally:
        driver.quit()
        print("Chrome fechado após processar CDMs.")
