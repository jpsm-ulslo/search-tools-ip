from _utils import iniciar_chrome, pesquisar_npdm_e_exportar
from _test import ler_cdms_excel, processar_cdms_individualmente
import os

def main():
    npdm = "Bombas de perfusão"
    print("Abrindo Chrome...")
    driver = iniciar_chrome()
    
    if False:
        try:
            # Etapa 1: Pesquisar NPDM e exportar resultados
            arquivo_resultado = pesquisar_npdm_e_exportar(driver, npdm)
            if arquivo_resultado:
                print(f"Arquivo de resultados salvo: {arquivo_resultado}")
            else:
                print("Falha ao exportar resultados.")
                return  # encerra se não exportou
        finally:
            driver.quit()
            print("Chrome fechado.")
    else:
        arquivo_resultado = "Referencias 29-08-2025 14_04.xlsx"

    # Etapa 2: Processar arquivo exportado
    if arquivo_resultado and os.path.exists(arquivo_resultado):
        try:
            cdms_info = ler_cdms_excel(arquivo_resultado)
            print(f"{len(cdms_info)} CDMs Comercializados encontrados.")
            processar_cdms_individualmente(cdms_info)
        except Exception as e:
            print(f"Erro durante o processamento dos CDMs: {e}")
        finally:
            print("Etapa de processamento finalizada.")
    else:
        print("Arquivo exportado não encontrado, etapa de processamento não será executada.")

if __name__ == "__main__":
    main()
