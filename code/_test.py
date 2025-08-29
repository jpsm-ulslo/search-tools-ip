import os
import time
import csv
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait, Select
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import TimeoutException, NoSuchElementException

from _utils import iniciar_chrome, digitar_simulado, aguardar_download

DOWNLOAD_DIR = os.getcwd()

def ler_cdms_csv(caminho_csv):
    """Ler arquivo CSV exportado pelo Infarmed e retornar lista de CDMs com Estado = Comercializado"""
    cdms = []
    with open(caminho_csv, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for linha in reader:
            estado = linha.get("Estado", "").strip()
            if estado.lower() == "comercializado":
                cdms.append({
                    "Referência": linha.get("Referência", "").strip(),
                    "Fabricante": linha.get("Fabricante", "").strip(),
                    "CDM": linha.get("CDM", "").strip(),
                    "Distribuidor": linha.get("Distribuidor", "").strip()
                })
    return cdms

def processar_cdms_individualmente(cdms_info):
    """Para cada CDM, acessar a página do Infarmed e buscar detalhes e folheto se disponível"""
    driver = iniciar_chrome(DOWNLOAD_DIR)
    wait = WebDriverWait(driver, 20)

    try:
        driver.get("https://www.infarmed.pt/web/infarmed/pesquisa-dispositivos")
        for idx, info in enumerate(cdms_info, 1):
            cdm = info["CDM"]
            print(f"[{idx}/{len(cdms_info)}] Processando CDM: {cdm}")
            try:
                # Campo de pesquisa individual
                cdm_input = wait.until(EC.presence_of_element_located(
                    (By.ID, "_SIDMPesquisaDispositivos_WAR_SIDMportlet_:pesquisaForm:infarmed5_input")
                ))
                cdm_input.clear()
                digitar_simulado(cdm_input, cdm)
                time.sleep(0.3)

                # Botão pesquisar
                search_btn = wait.until(EC.element_to_be_clickable(
                    (By.ID, "_SIDMPesquisaDispositivos_WAR_SIDMportlet_:pesquisaForm:pesquisarBtn")
                ))
                driver.execute_script("arguments[0].click();", search_btn)

                # Espera resultado da tabela individual
                linha_certa = wait.until(EC.presence_of_element_located(
                    (By.XPATH, f"//table[contains(@class,'ui-datatable')]//tr[td[3][text()='{cdm}']]")
                ))

                # Aqui poderia abrir detalhes e baixar folheto (similar ao que já estava em infarmed_utils)
                # info["Folheto"] = ...

            except TimeoutException:
                print(f"  - CDM {cdm} não encontrado ou sem resultados.")

    finally:
        driver.quit()
