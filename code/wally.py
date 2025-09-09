import os
import time
import shutil
import random
import re
import pandas as pd
from datetime import datetime
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import WebDriverWait, Select
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import TimeoutException, NoSuchElementException, WebDriverException

# ------------------------------
# Configurações
# ------------------------------
BASE_DIR = os.getcwd()
TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")
DOWNLOAD_DIR = os.path.join(BASE_DIR, f"resultados_{TIMESTAMP}")
os.makedirs(DOWNLOAD_DIR, exist_ok=True)

TIMEOUT_PESQUISA = 20
TIMEOUT_DOWNLOAD = 120
MAX_RETRIES_PDF = 3
MAX_RETRIES_SITE = 3

# ------------------------------
# Funções utilitárias
# ------------------------------
def iniciar_chrome(download_dir=DOWNLOAD_DIR):
    options = webdriver.ChromeOptions()
    prefs = {
        "download.default_directory": download_dir,
        "download.prompt_for_download": False,
        "directory_upgrade": True
    }
    options.add_experimental_option("prefs", prefs)
    driver = webdriver.Chrome(options=options)
    driver.set_page_load_timeout(30)
    return driver

def digitar_simulado(elemento, texto, delay_min=0.05, delay_max=0.15):
    for letra in texto:
        elemento.send_keys(letra)
        time.sleep(random.uniform(delay_min, delay_max))

def sanitize_filename(name):
    return re.sub(r'[\\/*?:"<>|]', "_", name)

def aguardar_download(diretorio, arquivos_antes, timeout=TIMEOUT_DOWNLOAD):
    inicio = time.time()
    while time.time() - inicio < timeout:
        arquivos_atuais = set(os.listdir(diretorio))
        novos = arquivos_atuais - arquivos_antes
        novos_validos = [f for f in novos if not f.endswith((".tmp", ".crdownload"))]
        if novos_validos:
            return novos_validos[0]
        time.sleep(0.5)
    return None

# ------------------------------
# Função de baixar folheto PDF
# ------------------------------
def baixar_folheto(driver, modal, cdm, max_tentativas=MAX_RETRIES_PDF):
    for tentativa in range(1, max_tentativas + 1):
        try:
            folheto_icon = modal.find_element(By.CSS_SELECTOR,
                "span.icon-download-alt[title='Folheto Informativo']")
            driver.execute_script("arguments[0].click();", folheto_icon)
            arquivos_antes = set(os.listdir(DOWNLOAD_DIR))
            arquivo_pdf = aguardar_download(DOWNLOAD_DIR, arquivos_antes)
            if arquivo_pdf:
                destino = os.path.join(DOWNLOAD_DIR, f"{sanitize_filename(cdm)}_folheto.pdf")
                shutil.move(os.path.join(DOWNLOAD_DIR, arquivo_pdf), destino)
                print(f"PDF baixado: {destino}")
                return True
        except (NoSuchElementException, TimeoutException):
            print(f"  - PDF não disponível na tentativa {tentativa} para CDM {cdm}")
        time.sleep(1)
    print(f"❌ Falha: não foi possível baixar PDF do CDM {cdm}")
    return False

# ------------------------------
# Função de pesquisa e exportação de NPDM
# ------------------------------
def pesquisar_npdm_e_exportar(driver, npdm="Bombas de perfusão",
                               pasta_destino=DOWNLOAD_DIR,
                               timeout_download=TIMEOUT_DOWNLOAD,
                               max_tentativas=MAX_RETRIES_SITE):
    for tentativa in range(1, max_tentativas + 1):
        print(f"[Tentativa {tentativa}/{max_tentativas}] Pesquisando NPDM: {npdm}")
        try:
            driver.get("https://www.infarmed.pt/web/infarmed/pesquisa-dispositivos")
            wait = WebDriverWait(driver, TIMEOUT_PESQUISA)

            npdm_input = wait.until(EC.presence_of_element_located(
                (By.ID, "_SIDMPesquisaDispositivos_WAR_SIDMportlet_:pesquisaForm:infarmed6_input")
            ))
            npdm_input.clear()
            digitar_simulado(npdm_input, npdm)
            npdm_input.send_keys("\n")

            # Tentar filtrar "Comercializado"
            try:
                estado_select = Select(wait.until(
                    EC.element_to_be_clickable(
                        (By.ID, "_SIDMPesquisaDispositivos_WAR_SIDMportlet_:pesquisaForm:estadoComercializacao")
                    )
                ))
                estado_select.select_by_visible_text("Comercializado")
            except TimeoutException:
                print("⚠️ Campo de estado não encontrado, continuando...")

            search_btn = wait.until(EC.element_to_be_clickable(
                (By.ID, "_SIDMPesquisaDispositivos_WAR_SIDMportlet_:pesquisaForm:pesquisarBtn")
            ))
            driver.execute_script("arguments[0].click();", search_btn)
            time.sleep(2)

            export_btn = wait.until(EC.element_to_be_clickable(
                (By.ID, "_SIDMPesquisaDispositivos_WAR_SIDMportlet_:pesquisaForm:exportarBtn")
            ))
            arquivos_antes = set(os.listdir(pasta_destino))
            driver.execute_script("arguments[0].click();", export_btn)

            arquivo_baixado = aguardar_download(pasta_destino, arquivos_antes, timeout=timeout_download)
            if not arquivo_baixado:
                print(f"❌ Download não detectado dentro de {timeout_download}s.")
                continue

            caminho_origem = os.path.join(pasta_destino, arquivo_baixado)
            if not os.path.exists(caminho_origem):
                print(f"❌ Arquivo baixado não encontrado: {caminho_origem}")
                continue

            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            nome_final = f"{sanitize_filename(npdm)}_{timestamp}.xlsx"
            destino_final = os.path.join(pasta_destino, nome_final)
            shutil.move(caminho_origem, destino_final)
            print(f"Download finalizado e renomeado: {destino_final}")
            return destino_final

        except (TimeoutException, WebDriverException) as e:
            print(f"Erro na tentativa {tentativa}: {e}")
            time.sleep(2)

    print(f"❌ Falha final: não foi possível pesquisar/exportar NPDM '{npdm}'")
    return None

# ------------------------------
# Processamento de CDMs
# ------------------------------
def processar_cdms(cdms_info, relatorio_path=None):
    driver = iniciar_chrome()
    wait = WebDriverWait(driver, TIMEOUT_PESQUISA)
    resultados = []

    try:
        driver.get("https://www.infarmed.pt/web/infarmed/pesquisa-dispositivos")
        for idx, info in enumerate(cdms_info, 1):
            cdm = info.get("CDM")
            print(f"[{idx}/{len(cdms_info)}] Processando CDM: {cdm}")
            resultado = info.copy()
            resultado["Status"] = "Não processado"
            resultado["Folheto_Baixado"] = False

            try:
                cdm_input = wait.until(EC.presence_of_element_located(
                    (By.ID, "_SIDMPesquisaDispositivos_WAR_SIDMportlet_:pesquisaForm:infarmed5_input")
                ))
                cdm_input.clear()
                digitar_simulado(cdm_input, cdm)

                search_btn = wait.until(EC.element_to_be_clickable(
                    (By.ID, "_SIDMPesquisaDispositivos_WAR_SIDMportlet_:pesquisaForm:pesquisarBtn")
                ))
                driver.execute_script("arguments[0].click();", search_btn)
                time.sleep(1)

                try:
                    detalhes_icon = wait.until(EC.element_to_be_clickable(
                        (By.CSS_SELECTOR, "div.icon-file-alt[title='Detalhes']")
                    ))
                    driver.execute_script("arguments[0].click();", detalhes_icon)
                    modal = wait.until(EC.visibility_of_element_located(
                        (By.ID, "_SIDMPesquisaDispositivos_WAR_SIDMportlet_:detalheForm:detalheDialog")
                    ))
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
                    except TimeoutException:
                        pass

                except TimeoutException:
                    print(f"  - CDM {cdm} não encontrado ou sem resultados.")
                    resultado["Status"] = "Não encontrado"

            except (TimeoutException, WebDriverException) as e:
                print(f"  - Erro ao pesquisar CDM {cdm}: {e}")
                resultado["Status"] = "Erro"

            resultados.append(resultado)

    finally:
        driver.quit()
        print("Chrome fechado após processar CDMs.")

    if relatorio_path is None:
        relatorio_path = os.path.join(DOWNLOAD_DIR, "relatorio_final.xlsx")
    df = pd.DataFrame(resultados)
    df.to_excel(relatorio_path, index=False)
    print(f"📄 Relatório salvo em: {relatorio_path}")
    return relatorio_path

# ------------------------------
# Função principal
# ------------------------------
def main():
    npdm = "Bicicletas para usos fisioterapêuticos e/ou diagnósticos"
    print("Abrindo Chrome...")
    driver = iniciar_chrome()
    arquivo_resultado = None

    try:
        arquivo_resultado = pesquisar_npdm_e_exportar(
            driver, npdm=npdm, pasta_destino=DOWNLOAD_DIR,
            timeout_download=TIMEOUT_DOWNLOAD, max_tentativas=MAX_RETRIES_SITE
        )
        if arquivo_resultado:
            print(f"Arquivo de resultados salvo: {arquivo_resultado}")
        else:
            print("Falha ao exportar resultados.")
            return
    finally:
        driver.quit()
        print("Chrome fechado após pesquisa.")

    if arquivo_resultado and os.path.exists(arquivo_resultado):
        try:
            df = pd.read_excel(arquivo_resultado)
            cdms_info = []
            for _, row in df.iterrows():
                cdms_info.append({
                    "CDM": str(row.get("CDM", "")).strip(),
                    "Referência": str(row.get("Referência", "")).strip(),
                    "Fabricante": str(row.get("Fabricante", "")).strip(),
                    "Distribuidor": str(row.get("Distribuidor", "")).strip(),
                })
            print(f"{len(cdms_info)} CDMs Comercializados encontrados.")
            processar_cdms(cdms_info)
        except Exception as e:
            print(f"Erro durante o processamento dos CDMs: {e}")
    else:
        print("Arquivo exportado não encontrado, processamento não será executado.")

# ------------------------------
# Entrypoint
# ------------------------------
if __name__ == "__main__":
    main()
