import javax.swing.JOptionPane;

public class CalculadoraSwing {

    private int resultat = 0;

    public void demanarISumar() {
        try {
            String aStr = JOptionPane.showInputDialog(null, "Introdueix el primer número:");
            String bStr = JOptionPane.showInputDialog(null, "Introdueix el segon número:");

            int a = Integer.parseInt(aStr);
            int b = Integer.parseInt(bStr);

            resultat = a + b;

            JOptionPane.showMessageDialog(null, "La suma és: " + resultat);
        } catch (Exception e) {
            JOptionPane.showMessageDialog(null, "Error en l'entrada: " + e.getMessage());
            resultat = 0;
        }
    }

    public int obtenirResultat() {
        return resultat;
    }
}
