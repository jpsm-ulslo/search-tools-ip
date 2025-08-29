import os
import time
import random
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import WebDriverWait, Select
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import TimeoutException

# Diretório onde o script principal roda
DOWNLOAD_DIR = os.getcwd()

def digitar_simulado(elemento, texto, delay_min=0.1, delay_max=0.3):
    for letra in texto:
        elemento.send_keys(letra)
        time.sleep(random.uniform(delay_min, delay_max))

def aguardar_download(diretorio, arquivos_antes, timeout=60):
    inicio = time.time()
    while time.time() - inicio < timeout:
        arquivos_atuais = set(os.listdir(diretorio))
        novos = arquivos_atuais - arquivos_antes
        novos_validos = [f for f in novos if not f.endswith(".tmp") and f != os.path.basename(__file__)]
        tmp_exist = any(f.endswith(".tmp") for f in arquivos_atuais)
        if novos_validos and not tmp_exist:
            return novos_validos[0]
        time.sleep(0.5)
    return None

def iniciar_chrome(download_dir=DOWNLOAD_DIR):
    options = webdriver.ChromeOptions()
    prefs = {
        "download.default_directory": download_dir,
        "download.prompt_for_download": False,
        "directory_upgrade": True
    }
    options.add_experimental_option("prefs", prefs)
    driver = webdriver.Chrome(options=options)
    return driver

def pesquisar_npdm_e_exportar(driver, npdm="Bombas de perfusão"):
    wait = WebDriverWait(driver, 20)
    driver.get("https://www.infarmed.pt/web/infarmed/pesquisa-dispositivos")
    print(f"Digitando NPDM: {npdm}")
    
    npdm_input = wait.until(EC.presence_of_element_located(
        (By.ID, "_SIDMPesquisaDispositivos_WAR_SIDMportlet_:pesquisaForm:infarmed6_input")
    ))
    npdm_input.clear()
    digitar_simulado(npdm_input, npdm)
    time.sleep(0.5)
    npdm_input.send_keys(Keys.ENTER)
    
    try:
        estado_select = Select(wait.until(
            EC.element_to_be_clickable(
                (By.ID, "_SIDMPesquisaDispositivos_WAR_SIDMportlet_:pesquisaForm:estadoComercializacao")
            )
        ))
        estado_select.select_by_visible_text("Comercializado")
    except TimeoutException:
        print("Campo de estado não encontrado, continuando...")

    search_btn = wait.until(EC.element_to_be_clickable(
        (By.ID, "_SIDMPesquisaDispositivos_WAR_SIDMportlet_:pesquisaForm:pesquisarBtn")
    ))
    print("Clicando no botão Pesquisar...")
    driver.execute_script("arguments[0].click();", search_btn)
    time.sleep(2)

    # Exportar resultados
    export_btn = wait.until(EC.element_to_be_clickable(
        (By.ID, "_SIDMPesquisaDispositivos_WAR_SIDMportlet_:pesquisaForm:exportarBtn")
    ))
    print("Exportando resultados...")
    arquivos_antes = set(os.listdir(DOWNLOAD_DIR))
    driver.execute_script("arguments[0].click();", export_btn)
    arquivo_baixado = aguardar_download(DOWNLOAD_DIR, arquivos_antes)
    
    if arquivo_baixado:
        print(f"Download finalizado: {os.path.join(DOWNLOAD_DIR, arquivo_baixado)}")
        return os.path.join(DOWNLOAD_DIR, arquivo_baixado)
    else:
        print("Download não detectado dentro do tempo limite.")
        return None
