import scvi 
import scanpy as sc
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("input")
parser.add_argument("obsoutput")
args = parser.parse_args()

adata = sc.read_10x_h5(args.input)
print(adata)
scvi.model.SCVI.setup_anndata(adata)
vae = scvi.model.SCVI(adata)
try:
    vae.train(early_stopping = True, accelerator = 'gpu', batch_size =127)
except:
    vae.train(early_stopping = True, accelerator = 'gpu', batch_size =126)
    
solo = scvi.external.SOLO.from_scvi_model(vae)
try:
    solo.train(batch_size =127)
except:
    solo.train(batch_size =126)

preds = solo.predict()

preds['prediction'] = solo.predict(soft = False)

preds['solo_dif'] = preds.doublet - preds.singlet

doublets = preds[(preds.prediction == 'doublet') & (preds.solo_dif > 1)]

adata.obs['solo_doublet'] = adata.obs.index.isin(doublets.index)

adata.obs['solo_score'] = preds['solo_dif']

adata.obs.to_csv(args.obsoutput)
